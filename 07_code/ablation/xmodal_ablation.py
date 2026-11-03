"""
Input-level ablation framework: what information in spatial RNA transfers to spatial protein?

One fixed, deliberately simple predictor (ridge regression, one model per protein) is trained
to predict each protein from RNA-derived inputs. Every experiment changes ONLY the inputs:
which genes the model sees, how much spatial context it sees, whether that context is real
or shuffled, and which tissue regions it was trained on. Because the model never changes,
any change in held-out accuracy is attributable to the information that was removed.

Experiments (all return tidy pandas DataFrames, one row per protein x condition):
  1. variance_ladder     nested inputs: cell type -> + cognate gene -> + all RNA -> + neighbour RNA
  2. gene_ablation       drop cognate genes / random genes of equal size / named gene groups
  3. radius_sweep        neighbourhood radius 0..R um
  4. context_controls    real neighbourhood vs neighbours shuffled within cell type vs RNA shuffled
  5. region_holdout      train on all regions except one, test on the held-out region
  6. residual_map        protein minus RNA-predicted protein, with a within-cell-type Moran's I

Splits are always by sample (patient/slide), never by cell, so neighbouring cells cannot
leak across train and test. If there is only one sample, spatial blocks are used instead.

Inputs are plain arrays so any AnnData / SpatialData object can be adapted in a few lines
(see `from_anndata`). Run `python run_simulation_demo.py` for a self-check on simulated data.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder, StandardScaler

ALPHAS = np.logspace(-2, 4, 13)


# --------------------------------------------------------------------------------------
# data container
# --------------------------------------------------------------------------------------
@dataclass
class PairedSpatial:
    rna: np.ndarray                 # cells x genes, log1p-normalised
    protein: np.ndarray             # cells x proteins, arcsinh/z-scored per sample
    coords: np.ndarray              # cells x 2, micrometres
    cell_type: np.ndarray           # cells, str labels
    sample: np.ndarray              # cells, patient / slide id (split unit)
    genes: list[str]
    proteins: list[str]
    cognate: dict[str, list[str]] = field(default_factory=dict)   # protein -> gene(s)
    region: np.ndarray | None = None                              # cells, region / niche label
    gene_groups: dict[str, list[str]] = field(default_factory=dict)

    def __post_init__(self):
        if issparse(self.rna):
            self.rna = self.rna.toarray()
        self.rna = np.asarray(self.rna, dtype=float)
        self.protein = np.asarray(self.protein, dtype=float)
        self.coords = np.asarray(self.coords, dtype=float)
        self.cell_type = np.asarray(self.cell_type).astype(str)
        self.sample = np.asarray(self.sample).astype(str)
        self.gidx = {g: i for i, g in enumerate(self.genes)}


def from_anndata(ad_rna, protein_key="protein", spatial_key="spatial", celltype_key="cell_type",
                 sample_key="sample", region_key=None, cognate=None):
    """Build PairedSpatial from one AnnData with RNA in X and protein in obsm[protein_key]
    (a DataFrame, as written by Xenium/SpatialData readers or muon)."""
    prot = ad_rna.obsm[protein_key]
    prot_names = list(prot.columns) if hasattr(prot, "columns") else [f"p{i}" for i in range(prot.shape[1])]
    return PairedSpatial(
        rna=ad_rna.X, protein=np.asarray(prot), coords=ad_rna.obsm[spatial_key],
        cell_type=ad_rna.obs[celltype_key].values, sample=ad_rna.obs[sample_key].values,
        genes=list(ad_rna.var_names), proteins=prot_names, cognate=cognate or {},
        region=None if region_key is None else ad_rna.obs[region_key].values)


# --------------------------------------------------------------------------------------
# feature builders (the only thing experiments change)
# --------------------------------------------------------------------------------------
def radius_graph(coords, sample, radius):
    """Row-normalised neighbour matrix (excluding self) within `radius` um, within sample."""
    n = len(coords)
    rows, cols = [], []
    for s in np.unique(sample):
        idx = np.where(sample == s)[0]
        nn = NearestNeighbors(radius=radius).fit(coords[idx])
        _, nbrs = nn.radius_neighbors(coords[idx])
        for i, nb in zip(idx, nbrs):
            nb = idx[nb]
            nb = nb[nb != i]
            rows.extend([i] * len(nb))
            cols.extend(nb)
    W = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    deg = np.asarray(W.sum(1)).ravel()
    deg[deg == 0] = 1
    return csr_matrix(W.multiply(1 / deg[:, None]))


def shuffle_within(groups_a, groups_b, rng):
    """Permutation of cell indices that only swaps cells sharing (group_a, group_b)."""
    perm = np.arange(len(groups_a))
    key = pd.Series(list(zip(groups_a, groups_b)))
    for _, idx in key.groupby(key).groups.items():
        idx = np.asarray(idx)
        perm[idx] = rng.permutation(idx)
    return perm


class Features:
    """Builds design matrices from a PairedSpatial object under a given ablation."""

    def __init__(self, d: PairedSpatial, n_pcs=30, interact=True, seed=0):
        # interact=True lets cognate-gene and neighbourhood effects differ by the cell's own type
        # (e.g. IFNG-rich neighbours induce PD-L1 in tumour cells only). Still a linear model.
        self.d, self.n_pcs, self.interact, self.rng = d, n_pcs, interact, np.random.default_rng(seed)
        self.ct = OneHotEncoder(sparse_output=False).fit_transform(d.cell_type[:, None])
        self._W = {}
        self._red = {}   # frozen scaler+PCA per gene mask (unsupervised, fit once on intact RNA)

    def _reduce(self, X, mask_key):
        """Scale + PCA with transforms fitted ONCE per gene mask and then reused, so that
        permuted / shuffled inputs are projected on the same axes the model was trained on."""
        k = min(self.n_pcs, X.shape[1])
        if X.shape[1] <= k:
            return X
        if mask_key not in self._red:
            sc = StandardScaler().fit(X)
            self._red[mask_key] = (sc, PCA(k, random_state=0).fit(sc.transform(X)))
        sc, pca = self._red[mask_key]
        return pca.transform(sc.transform(X))

    def W(self, radius):
        if radius not in self._W:
            self._W[radius] = radius_graph(self.d.coords, self.d.sample, radius)
        return self._W[radius]

    def _ix(self, B):
        """Append (block x own-cell-type indicator) interactions."""
        if not self.interact:
            return B
        return np.hstack([B] + [B * self.ct[:, [k]] for k in range(self.ct.shape[1])])

    def _mask(self, drop_genes=(), keep_genes=None):
        mask = np.ones(self.d.rna.shape[1], bool)
        if keep_genes is not None:
            mask[:] = False
            mask[[self.d.gidx[g] for g in keep_genes if g in self.d.gidx]] = True
        for g in drop_genes:
            if g in self.d.gidx:
                mask[self.d.gidx[g]] = False
        return mask

    def rna_view(self, drop_genes=(), keep_genes=None, rna_perm=None):
        X = self.d.rna
        if rna_perm is not None:
            X = X[rna_perm]
        mask = np.ones(X.shape[1], bool)
        if keep_genes is not None:
            mask[:] = False
            mask[[self.d.gidx[g] for g in keep_genes if g in self.d.gidx]] = True
        for g in drop_genes:
            if g in self.d.gidx:
                mask[self.d.gidx[g]] = False
        return X[:, mask]

    def build(self, *, celltype=True, cognate_for=None, all_rna=False, nbr_radius=0,
              drop_genes=(), keep_genes=None, rna_perm=None, nbr_perm=None, nbr_composition=False):
        """Return feature matrix. `cognate_for` = protein name -> its cognate gene(s) only.
        `nbr_perm` permutes which cell's RNA a neighbour slot reads (context shuffle)."""
        blocks = []
        if celltype:
            blocks.append(self.ct)
        X = self.rna_view(drop_genes, keep_genes, rna_perm)
        if cognate_for is not None:
            cg = [g for g in self.d.cognate.get(cognate_for, []) if g in self.d.gidx and g not in drop_genes]
            if cg:
                src = self.d.rna if rna_perm is None else self.d.rna[rna_perm]
                blocks.append(self._ix(src[:, [self.d.gidx[g] for g in cg]]))
        key = self._mask(drop_genes, keep_genes).tobytes()
        if (all_rna or nbr_radius) and key not in self._red and X.shape[1] > 0:
            self._reduce(self.rna_view(drop_genes, keep_genes), key)   # fit on intact, unpermuted RNA
        if all_rna and X.shape[1] > 0:
            blocks.append(self._reduce(X, key))
        if nbr_radius and nbr_radius > 0 and X.shape[1] > 0:
            Xn = X if nbr_perm is None else X[nbr_perm]
            Z = self._reduce(Xn, key)
            blocks.append(self._ix(self.W(nbr_radius) @ Z))
            if nbr_composition:
                ct = self.ct if nbr_perm is None else self.ct[nbr_perm]
                blocks.append(self._ix(self.W(nbr_radius) @ ct))
        return np.hstack(blocks) if blocks else np.zeros((len(self.d.sample), 1))


# --------------------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------------------
def splits(d: PairedSpatial, n_blocks=4):
    """Leave-one-sample-out; if a single sample, leave-one-spatial-block-out (vertical strips)."""
    samples = np.unique(d.sample)
    if len(samples) > 1:
        for s in samples:
            yield s, d.sample != s, d.sample == s
    else:
        x = d.coords[:, 0]
        edges = np.quantile(x, np.linspace(0, 1, n_blocks + 1))
        blk = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, n_blocks - 1)
        for b in range(n_blocks):
            yield f"block{b}", blk != b, blk == b


def r2(y, yhat):
    ss = ((y - y.mean()) ** 2).sum()
    return np.nan if ss == 0 else 1 - ((y - yhat) ** 2).sum() / ss


def evaluate(d, make_X, label, folds=None, return_pred=False):
    """make_X(protein_name) -> feature matrix. Returns per-protein held-out R2 (mean over folds)."""
    folds = folds or list(splits(d))
    rows, pred = [], np.zeros_like(d.protein)
    for j, p in enumerate(d.proteins):
        X = make_X(p)
        y = d.protein[:, j]
        scores, cors = [], []
        for fold, tr, te in folds:
            sc = StandardScaler().fit(X[tr])
            m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[tr]), y[tr])
            yh = m.predict(sc.transform(X[te]))
            pred[te, j] = yh
            scores.append(r2(y[te], yh))
            cors.append(np.corrcoef(y[te], yh)[0, 1] if yh.std() > 0 and y[te].std() > 0 else np.nan)
        rows.append({"protein": p, "condition": label, "r2": float(np.nanmean(scores)),
                     "r2_sd": float(np.nanstd(scores)), "pearson": float(np.nanmean(cors))})
    out = pd.DataFrame(rows)
    return (out, pred) if return_pred else out


# --------------------------------------------------------------------------------------
# experiments
# --------------------------------------------------------------------------------------
def variance_ladder(d, F: Features, radius=30):
    """How much protein variance each information source adds, in a fixed order."""
    steps = [
        ("1_celltype", lambda p: F.build()),
        ("2_+cognate_gene", lambda p: F.build(cognate_for=p)),
        ("3_+all_rna", lambda p: F.build(cognate_for=p, all_rna=True)),
        (f"4_+nbr_rna_{radius}um", lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius)),
    ]
    res = pd.concat([evaluate(d, f, lab) for lab, f in steps])
    res["gain_vs_prev"] = res.groupby("protein")["r2"].diff()
    return res


def gene_ablation(d, F: Features, radius=30, n_random=5, seed=0):
    """Drop gene sets from ALL inputs (own + neighbour RNA) and measure the R2 loss."""
    rng = np.random.default_rng(seed)
    full = evaluate(d, lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius), "full")
    out = [full]
    cog_all = sorted({g for gs in d.cognate.values() for g in gs if g in d.gidx})
    out.append(evaluate(d, lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius,
                                             drop_genes=d.cognate.get(p, [])), "drop_own_cognate"))
    out.append(evaluate(d, lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius,
                                             drop_genes=cog_all), "drop_all_cognates"))
    non_cog = [g for g in d.genes if g not in set(cog_all)]
    for r in range(n_random):
        rand = list(rng.choice(non_cog, size=min(len(cog_all), len(non_cog)), replace=False))
        out.append(evaluate(d, lambda p, rand=rand: F.build(cognate_for=p, all_rna=True, nbr_radius=radius,
                                                            drop_genes=rand), f"drop_random_{r}"))
    for name, gs in d.gene_groups.items():
        out.append(evaluate(d, lambda p, gs=gs: F.build(cognate_for=p, all_rna=True, nbr_radius=radius,
                                                        drop_genes=gs), f"drop_group_{name}"))
    res = pd.concat(out)
    base = full.set_index("protein")["r2"]
    res["delta_vs_full"] = res["r2"] - res["protein"].map(base)
    return res


def radius_sweep(d, F: Features, radii=(0, 10, 20, 30, 50, 100)):
    return pd.concat([evaluate(d, lambda p, r=r: F.build(cognate_for=p, all_rna=True, nbr_radius=r),
                               f"radius_{r:03d}") for r in radii])


def context_controls(d, F: Features, radius=30, seed=0):
    """Is the neighbourhood gain about WHERE cells are, or only WHICH cell types are around?"""
    rng = np.random.default_rng(seed)
    own = lambda p: F.build(cognate_for=p, all_rna=True)
    real = lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius)
    # neighbour slots read RNA from a random cell of the SAME type in the same sample:
    # keeps local composition, destroys cell-state arrangement
    perm_ct = shuffle_within(d.sample, d.cell_type, rng)
    ct_shuf = lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius, nbr_perm=perm_ct)
    # neighbour slots read RNA from any random cell in the sample: destroys composition too
    perm_any = shuffle_within(d.sample, np.zeros(len(d.sample)), rng)
    any_shuf = lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius, nbr_perm=perm_any)
    # negative control: every RNA input comes from a random other cell (like arm G in H2)
    rna_shuf = lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius, rna_perm=perm_any)
    return pd.concat([evaluate(d, f, lab) for lab, f in [
        ("a_own_rna_only", own), ("b_+real_neighbourhood", real),
        ("c_+nbr_shuffled_within_celltype", ct_shuf), ("d_+nbr_shuffled_any", any_shuf),
        ("e_all_rna_shuffled", rna_shuf)]])


def region_holdout(d, F: Features, radius=30):
    """Spatial OOD transfer: train on all regions but one (all samples), test on that region."""
    if d.region is None:
        raise ValueError("region labels required")
    folds = [(f"holdout_{r}", d.region != r, d.region == r) for r in np.unique(d.region)]
    res = []
    for fold in folds:
        for lab, f in [("own_rna", lambda p: F.build(cognate_for=p, all_rna=True)),
                       ("own+nbr_rna", lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius))]:
            r = evaluate(d, f, lab, folds=[fold])
            r["heldout_region"] = fold[0].replace("holdout_", "")
            res.append(r)
    return pd.concat(res)


def morans_i_within_type(values, d: PairedSpatial, radius=30, n_perm=199, seed=0):
    """Moran's I of `values` computed only on edges between cells of the same type (per sample),
    with a within-type permutation p-value. Dependency-free test of within-type spatial structure (used for residual maps)."""
    rng = np.random.default_rng(seed)
    W = radius_graph(d.coords, d.sample, radius).tolil()
    rows, cols = W.nonzero()
    keep = d.cell_type[rows] == d.cell_type[cols]
    rows, cols = rows[keep], cols[keep]
    out = {}
    for ct in np.unique(d.cell_type):
        m = (d.cell_type[rows] == ct)
        r_, c_ = rows[m], cols[m]
        idx = np.where(d.cell_type == ct)[0]
        if len(r_) < 20:
            continue
        def I(v):
            z = v - v[idx].mean()
            return (z[r_] * z[c_]).sum() / len(r_) / (z[idx] ** 2).mean()
        obs = I(values)
        null = []
        for _ in range(n_perm):
            v = values.copy()
            v[idx] = rng.permutation(values[idx])
            null.append(I(v))
        out[ct] = (obs, (1 + np.sum(np.array(null) >= obs)) / (n_perm + 1))
    return out


def jitter_index(coords, sample, delta_um, rng, cell_type=None):
    """For each cell, the index of a random cell within delta_um (same sample; optionally same type).
    Simulates pairing/registration error: the cell's RNA is read from a nearby, wrong cell."""
    idx_out = np.arange(len(coords))
    if delta_um <= 0:
        return idx_out
    for s in np.unique(sample):
        idx = np.where(sample == s)[0]
        nn = NearestNeighbors(radius=delta_um).fit(coords[idx])
        _, nbrs = nn.radius_neighbors(coords[idx])
        for i, nb in zip(idx, nbrs):
            cand = idx[nb]
            cand = cand[cand != i]
            if cell_type is not None:
                cand = cand[cell_type[cand] == cell_type[i]]
            if len(cand):
                idx_out[i] = rng.choice(cand)
    return idx_out


def pairing_degradation(d, F: Features, deltas=(0, 10, 20, 50, 100), radius=30, same_type=False, seed=0):
    """E8: RNA of each cell replaced by RNA of a random cell within delta um (all features rebuilt)."""
    rng = np.random.default_rng(seed)
    out = []
    for de in deltas:
        perm = jitter_index(d.coords, d.sample, de, rng, d.cell_type if same_type else None)
        out.append(evaluate(d, lambda p, perm=perm: F.build(cognate_for=p, all_rna=True, nbr_radius=radius,
                                                            rna_perm=perm), f"jitter_{de:03d}um"))
    return pd.concat(out)


def residual_map(d, F: Features, radius=30):
    """Protein not explained by RNA (own + neighbourhood). Spatial structure left in the residual
    = candidate post-transcriptional / microenvironment-driven protein regulation."""
    res, pred = evaluate(d, lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=radius),
                         "full", return_pred=True)
    resid = d.protein - pred
    rows = []
    for j, p in enumerate(d.proteins):
        for ct, (I_obs, pval) in morans_i_within_type(resid[:, j], d, radius).items():
            I_raw, p_raw = morans_i_within_type(d.protein[:, j], d, radius, n_perm=99).get(ct, (np.nan, np.nan))
            rows.append({"protein": p, "cell_type": ct, "moran_I_protein": I_raw,
                         "moran_I_residual": I_obs, "p_residual": pval})
    return pd.DataFrame(rows), resid
