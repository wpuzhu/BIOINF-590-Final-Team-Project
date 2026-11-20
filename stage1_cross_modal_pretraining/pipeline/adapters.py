"""
Model adapters for Stage 1: one interface, so the same arms (C1–C4, A1–A3) run on every model.

  adapter = get_adapter(name, cfg)            # "pca", "nicheformer", "spatium"
  ids, mask = adapter.encode(values, feature_names, modality, arm, q, rng)
  Z = adapter.embed(ids, mask, arm)                       # frozen mode
  clf = adapter.finetune(ids, mask, y, train_idx, arm)    # fine-tune mode

Cross-modal input is the crux of Stage 1:
  Nicheformer (gene-token vocabulary) reading PROTEIN: each marker is written at the token of
    its encoding gene (protein_cognate_gene); markers without a gene (PTMs, panCK) get no
    token. Values are ranked against the dataset's per-marker mean, as Nicheformer ranks genes
    against the technology mean.
  Spatium (protein-token vocabulary) reading RNA: each gene with a cognate marker in Spatium's
    vocabulary is written at that marker's token; all other genes are dropped (e.g. 135 of
    5,001 Xenium genes).
Arms change only the mapping or the weights:
  pretrained · scratch (same architecture, random weights) · A1_reinit (re-initialise the
  embedding rows of the mapped tokens) · A2_random_map (map to random wrong tokens of the same
  vocabulary) · A3_overlap (keep a fraction q of correct mappings, put the rest on fresh tokens)

Model-specific calls are marked VERIFY: they follow the public docs as of Nov 2026
(Nicheformer: Hugging Face "theislab/Nicheformer" with trust_remote_code, tokenizer with a
technology mean, model.get_embeddings; Spatium: spatium.readthedocs.io FineTuneTokenizer and
zarr fine-tuning). Check them against the installed versions before the first real run.
"""
from __future__ import annotations

import numpy as np

ARMS = ("pretrained", "scratch", "A1_reinit", "A2_random_map", "A3_overlap")


class PCAAdapter:
    """Reference bar (Justin's benchmark): per-patient z-score -> PCA -> logistic regression."""
    name = "pca"

    def __init__(self, cfg):
        self.k = cfg.get("n_pcs", 50)

    def encode(self, values, names, modality, arm, q, rng, patient=None):
        Z = np.asarray(values, float).copy()
        if modality == "rna":
            Z = np.log1p(Z / (Z.sum(1, keepdims=True) + 1e-9) * 1e4)
        if patient is not None:
            for p in np.unique(patient):
                m = patient == p
                Z[m] = (Z[m] - Z[m].mean(0)) / (Z[m].std(0) + 1e-9)
        return Z, None

    def embed(self, X, mask, arm, train_idx=None):
        from sklearn.decomposition import PCA
        k = min(self.k, X.shape[1] - 1)
        fit_on = X if train_idx is None else X[train_idx]
        return PCA(k, random_state=0).fit(fit_on).transform(X)

    def finetune(self, *a, **k):
        raise NotImplementedError("PCA is frozen-only")


class _TokenModel:
    """Shared mapping logic for gene- or protein-token models."""
    vocab: dict

    def map_features(self, names, modality, arm, q, rng):
        """Vocabulary id per input feature (-1 = dropped) and the ids whose rows A1 re-initialises."""
        tok = np.array([self.vocab.get(self.translate(n, modality), -1) for n in names])
        mapped = np.where(tok >= 0)[0]
        if arm == "A2_random_map" and len(mapped):
            pool = np.setdiff1d(np.array(list(self.vocab.values())), tok[mapped])
            tok[mapped] = rng.choice(pool, size=len(mapped), replace=False)
        elif arm == "A3_overlap" and len(mapped):
            drop = rng.choice(mapped, size=int(round((1 - q) * len(mapped))), replace=False)
            spare = self.spare_tokens(len(drop))
            tok[drop] = spare
        reinit = tok[mapped] if arm == "A1_reinit" else np.array([], int)
        return tok, reinit


class NicheformerAdapter(_TokenModel):
    name = "nicheformer"

    def __init__(self, cfg):
        from transformers import AutoModelForMaskedLM, AutoTokenizer  # VERIFY versions
        self.cfg = cfg
        self.hf_id = cfg.get("hf_id", "theislab/Nicheformer")
        self.tokenizer = AutoTokenizer.from_pretrained(self.hf_id, trust_remote_code=True)
        self.vocab = self.tokenizer.get_vocab()
        self.model = AutoModelForMaskedLM.from_pretrained(self.hf_id, trust_remote_code=True)
        self.max_len = cfg.get("max_len", 1500)
        self.cognate = {}

    def translate(self, name, modality):
        # protein -> encoding gene symbol (from the harmonised protein_cognate_gene column)
        return self.cognate.get(name, "") if modality == "protein" else name

    def spare_tokens(self, n):
        # genes the target never expresses = tokens never used in this dataset (fresh meaning)
        used = set(self._used)
        free = [i for g, i in self.vocab.items() if i not in used and not g.startswith("[")]
        return np.array(free[:n])

    def encode(self, values, names, modality, arm, q, rng, patient=None, cognate=None):
        if cognate is not None:
            self.cognate = dict(zip(names, cognate))
        self._used = [self.vocab.get(self.translate(n, modality), -1) for n in names]
        tok, self._reinit = self.map_features(names, modality, arm, q, rng)
        V = np.asarray(values, float)
        V = np.expm1(V) if modality == "protein" else V            # undo arcsinh-like compression
        rel = V / (V.mean(0) + 1e-9)                                 # rank against dataset mean
        keep = tok >= 0
        rel, tok = rel[:, keep], tok[keep]
        order = np.argsort(-rel, 1)[:, : self.max_len]
        ids = tok[order]
        mask = np.take_along_axis(rel, order, 1) > 0
        # VERIFY: Nicheformer prepends context tokens (modality / organism / assay); take them
        # from the tokenizer for 'spatial', 'human' and the closest assay (e.g. Xenium).
        ctx = np.array(self.cfg.get("context_token_ids", []), int)
        if len(ctx):
            ids = np.c_[np.tile(ctx, (len(ids), 1)), ids]
            mask = np.c_[np.ones((len(ids), len(ctx)), bool), mask]
        return ids, mask

    def _weights(self, arm):
        import copy

        import torch
        m = copy.deepcopy(self.model)
        if arm == "scratch":
            m = type(self.model)(self.model.config)                 # same architecture, random init
        elif arm == "A1_reinit" and len(self._reinit):
            emb = m.get_input_embeddings() if hasattr(m, "get_input_embeddings") else m.embeddings  # VERIFY
            with torch.no_grad():
                torch.nn.init.normal_(emb.weight[torch.tensor(self._reinit)], std=0.02)
        return m

    def embed(self, ids, mask, arm, batch=256, device="cuda"):
        import torch
        m = self._weights(arm).to(device).eval()
        out = []
        with torch.no_grad():
            for i in range(0, len(ids), batch):
                x = torch.tensor(ids[i:i + batch], device=device)
                am = torch.tensor(mask[i:i + batch], device=device)
                h = m.get_embeddings(x, attention_mask=am, layer=-1)  # VERIFY: (B, L, d) or (B, d)
                if h.dim() == 3:
                    h = (h * am[..., None]).sum(1) / am.sum(1, keepdim=True).clamp(min=1)
                out.append(h.float().cpu().numpy())
        return np.vstack(out)

    def finetune(self, ids, mask, y, train_idx, arm, n_classes, epochs=10, lr=1e-4, batch=128, device="cuda"):
        import torch
        import torch.nn as nn
        m = self._weights(arm).to(device)
        d = m.config.dim_model if hasattr(m.config, "dim_model") else m.config.hidden_size  # VERIFY
        head = nn.Linear(d, n_classes).to(device)
        opt = torch.optim.AdamW(list(m.parameters()) + list(head.parameters()), lr=lr, weight_decay=0.01)

        def forward(i):
            x = torch.tensor(ids[i], device=device)
            am = torch.tensor(mask[i], device=device)
            h = m.get_embeddings(x, attention_mask=am, layer=-1)
            if h.dim() == 3:
                h = (h * am[..., None]).sum(1) / am.sum(1, keepdim=True).clamp(min=1)
            return head(h)
        yt = torch.tensor(y, device=device)
        rng = np.random.default_rng(0)
        for _ in range(epochs):
            m.train()
            for b in np.array_split(rng.permutation(train_idx), max(1, len(train_idx) // batch)):
                loss = nn.functional.cross_entropy(forward(b), yt[b])
                opt.zero_grad(); loss.backward(); opt.step()
        m.eval()

        def predict(idx):
            with torch.no_grad():
                return np.concatenate([forward(b).argmax(1).cpu().numpy() for b in np.array_split(idx, max(1, len(idx) // 512))])
        return predict


class SpatiumAdapter(_TokenModel):
    """Spatium (spatium.readthedocs.io). Its loaders read zarr; we write a zarr per (dataset,
    arm) with protein-token ids and values and call its encoder for embeddings."""
    name = "spatium"

    def __init__(self, cfg):
        import spatium  # noqa: F401  VERIFY package / import name (repo: github.com/ploughhh/Spatium)
        self.cfg = cfg
        self.vocab = dict(np.load(cfg["protein_vocab"], allow_pickle=True).item())  # protein_id_map, VERIFY
        self.gene2marker = {}

    def translate(self, name, modality):
        return self.gene2marker.get(name, "") if modality == "rna" else name

    def spare_tokens(self, n):
        used = set(self._used)
        return np.array([i for i in self.vocab.values() if i not in used][:n])

    def encode(self, values, names, modality, arm, q, rng, patient=None, cognate=None, marker_gene=None):
        if marker_gene is not None:                      # {marker: gene} -> {gene: marker}
            self.gene2marker = {g: mk for mk, g in marker_gene.items() if g}
        self._used = [self.vocab.get(self.translate(n, modality), -1) for n in names]
        tok, self._reinit = self.map_features(names, modality, arm, q, rng)
        keep = tok >= 0
        V = np.asarray(values, float)[:, keep]
        if modality == "rna":
            V = np.log1p(V / (np.asarray(values, float).sum(1, keepdims=True) + 1e-9) * 1e4)
        print(f"  spatium reads {keep.sum()} of {len(names)} {modality} features")
        return np.tile(tok[keep], (len(V), 1)), V

    def embed(self, ids, values, arm, device="cuda"):
        raise NotImplementedError(
            "VERIFY: load the pretrained spaProFormer (or a random-init copy for 'scratch'; re-init "
            "rows self._reinit for A1), feed (ids, values) through its FineTuneTokenizer-compatible "
            "zarr, and return the cell embedding. See spatium.readthedocs.io 'fine_tune_train'.")

    def finetune(self, *a, **k):
        raise NotImplementedError("VERIFY: use Spatium's cell_type_prediction fine-tuning with the "
                                  "zarr written by encode(); split_by_file=True with train/val zarr.")


def get_adapter(name, cfg):
    return {"pca": PCAAdapter, "nicheformer": NicheformerAdapter, "spatium": SpatiumAdapter}[name](cfg)
