"""
Paired spatial RNA + protein simulator for Stage 2 (and the Stage 1 sanity checks).

Every cell has a location, a niche, a coarse type, a subtype (the fine label to predict),
an RNA profile and a protein profile. Four scenarios set WHERE the subtype information lives,
which is exactly the informative vs transferable distinction:

  null                  subtype is visible in protein only; RNA carries no subtype signal,
                        so RNA adds nothing (teacher = protein-only model).
  informative_only      subtype is visible ONLY in RNA (RNA-specific). A model that sees RNA
                        gains (informative) but nothing can be taught to a protein-only student
                        (not transferable).
  transferable_protein  subtype is visible strongly in RNA and weakly, spread across many
                        markers, in protein. Few labels cannot find the weak protein direction;
                        an RNA teacher labelling many unlabelled cells can (transferable).
  transferable_spatial  subtype is invisible in the cell's own protein but tied to the niche;
                        the niche is readable from the protein of neighbours. RNA reveals the
                        subtype directly; the student can only learn it through space.

Units: micrometres. Default 1,500 cells per 1 mm x 1 mm slide, one slide per patient.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.neighbors import NearestNeighbors

SCENARIOS = ("null", "informative_only", "transferable_protein", "transferable_spatial")


@dataclass
class SimConfig:
    scenario: str = "transferable_protein"
    n_patients: int = 10
    cells_per_patient: int = 1500
    side_um: float = 1000.0
    n_niches: int = 4
    n_coarse: int = 4             # coarse types; each has 2 subtypes -> 8 fine labels
    n_genes: int = 150
    n_markers: int = 20
    niche_scale_um: float = 180.0  # spatial length scale of niches
    niche_pref: float = 0.55       # P(dominant coarse type | niche)
    rna_subtype: float = 0.25      # subtype effect size in RNA (per gene sd; 150 genes)
    prot_subtype: dict = field(default_factory=lambda: {
        "null": 0.8, "informative_only": 0.0, "transferable_protein": 0.45, "transferable_spatial": 0.0})
    spatial_subtype_purity: float = 0.9
    batch_sd: float = 0.3
    seed: int = 0


@dataclass
class PairedData:
    coords: np.ndarray      # (n, 2)
    patient: np.ndarray     # (n,)
    niche: np.ndarray
    coarse: np.ndarray
    fine: np.ndarray        # label to predict: coarse*2 + subtype
    rna: np.ndarray         # (n, n_genes) log-scale expression
    protein: np.ndarray     # (n, n_markers)
    cfg: SimConfig


def _smooth_fields(rng, xy, k, scale, n_freq=12):
    """k independent smooth random fields (random Fourier features)."""
    w = rng.normal(0, 1 / scale, size=(k, n_freq, 2))
    ph = rng.uniform(0, 2 * np.pi, size=(k, n_freq))
    return np.stack([np.cos(xy @ w[j].T + ph[j]).sum(1) for j in range(k)], 1)


def simulate(cfg: SimConfig) -> PairedData:
    assert cfg.scenario in SCENARIOS, cfg.scenario
    rng = np.random.default_rng(cfg.seed)
    C, G, P = cfg.n_coarse, cfg.n_genes, cfg.n_markers
    # shared biology (same across patients)
    M_rna = rng.normal(0, 1.0, (C, G))
    S_rna = rng.normal(0, cfg.rna_subtype, (C, G)) * (0.0 if cfg.scenario == "null" else 1.0)
    M_pro = rng.normal(0, 1.5, (C, P))
    S_pro = rng.normal(0, 1.0, (C, P))
    S_pro *= cfg.prot_subtype[cfg.scenario] / np.sqrt((S_pro ** 2).mean(1, keepdims=True))
    niche_type = np.full((cfg.n_niches, C), (1 - cfg.niche_pref) / (C - 1))
    for k in range(cfg.n_niches):
        niche_type[k, k % C] = cfg.niche_pref
    # niches -> subtype rule for the spatial scenario (niche parity)
    niche_sub = np.arange(cfg.n_niches) % 2

    out = {k: [] for k in ("coords", "patient", "niche", "coarse", "fine", "rna", "protein")}
    for p in range(cfg.n_patients):
        n = cfg.cells_per_patient
        xy = rng.uniform(0, cfg.side_um, (n, 2))
        niche = _smooth_fields(rng, xy, cfg.n_niches, cfg.niche_scale_um).argmax(1)
        coarse = np.array([rng.choice(C, p=niche_type[k]) for k in niche])
        if cfg.scenario == "transferable_spatial":
            flip = rng.random(n) > cfg.spatial_subtype_purity
            sub = np.where(flip, 1 - niche_sub[niche], niche_sub[niche])
        else:
            sub = rng.integers(0, 2, n)
        sgn = (2 * sub - 1)[:, None] / 2.0
        b_rna = rng.normal(0, cfg.batch_sd, G)
        b_pro = rng.normal(0, cfg.batch_sd, P)
        rna = M_rna[coarse] + sgn * S_rna[coarse] + b_rna + rng.normal(0, 1, (n, G))
        pro = M_pro[coarse] + sgn * S_pro[coarse] + b_pro + rng.normal(0, 1, (n, P))
        for k, v in zip(out, (xy, np.full(n, p), niche, coarse, coarse * 2 + sub, rna, pro)):
            out[k].append(v)
    cat = {k: np.concatenate(v) for k, v in out.items()}
    return PairedData(cfg=cfg, **cat)


# --------------------------------------------------------------------------------------
def neighbour_mean(X, coords, patient, radius=40.0, include_self=False):
    """Mean of X over neighbours within `radius` um on the same slide (the 1-hop graph message)."""
    out = np.zeros_like(X)
    for p in np.unique(patient):
        idx = np.where(patient == p)[0]
        nn = NearestNeighbors(radius=radius).fit(coords[idx])
        nb = nn.radius_neighbors(coords[idx], return_distance=False)
        for i, js in enumerate(nb):
            if not include_self:
                js = js[js != i]
            out[idx[i]] = X[idx[js]].mean(0) if len(js) else X[idx[i]]
    return out


def jitter_pairing(coords, patient, delta, seed=0):
    """Index map that pairs each cell with the cell nearest to a point `delta` um away on the
    same slide: a simple model of adjacent-section (registration-based) pairing."""
    rng = np.random.default_rng(seed)
    m = np.arange(len(coords))
    if delta <= 0:
        return m
    for p in np.unique(patient):
        idx = np.where(patient == p)[0]
        ang = rng.uniform(0, 2 * np.pi, len(idx))
        tgt = coords[idx] + delta * np.c_[np.cos(ang), np.sin(ang)]
        nn = NearestNeighbors(n_neighbors=1).fit(coords[idx])
        m[idx] = idx[nn.kneighbors(tgt, return_distance=False)[:, 0]]
    return m


def graded_null(X, mode, coarse=None, niche=None, patient=None, seed=0):
    """Teacher-input nulls. 'all' shuffle across cells; 'within_type' shuffle within coarse type
    (destroys subtype info, keeps type info); 'within_niche' shuffle within niche (keeps niche
    info only); 'type_average' replace by the coarse-type mean. All within patient."""
    rng = np.random.default_rng(seed)
    Y = X.copy()
    if mode == "none":
        return Y
    groups = {"all": np.zeros(len(X), int), "within_type": coarse, "within_niche": niche,
              "type_average": coarse}[mode]
    for p in np.unique(patient):
        for g in np.unique(groups):
            idx = np.where((patient == p) & (groups == g))[0]
            if mode == "type_average":
                Y[idx] = X[idx].mean(0)
            else:
                Y[idx] = X[rng.permutation(idx)]
    return Y
