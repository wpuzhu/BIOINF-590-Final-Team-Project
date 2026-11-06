"""
Simulated biology for the H1 transfer-ablation testbed.

One latent cell-state model generates BOTH modalities:

    cell type  ->  latent programmes z (lineage programmes strong, subtype programmes weak)
    RNA        ->  G genes = loadings @ z + noise, Poisson-sampled with library size (sparse counts)
    protein    ->  P markers; marker j has a cognate gene g(j). Its level is
                   rho_j * (cognate gene's noise-free expression)            <- shared via identity
                 + (1 - rho_j) * (protein-specific readout of the same z)    <- shared via structure only
                 + measurement noise + per-patient batch scale

rho_j is the RNA-protein coupling of marker j. Scenarios set it differently, so the ground
truth about WHAT can transfer from an RNA-pretrained model to protein is known:

  coupled     rho = 0.9 for all markers            identity should transfer
  decoupled   rho = 0.15 for all markers           identity is misleading; only latent structure is shared
  realistic   lineage markers rho = 0.85, state markers rho = 0.3,
              and 4 of 16 markers have NO gene partner (DNA stain / pan-keratin-like)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

N_LINEAGE, N_SUB = 4, 3                      # 12 fine cell types
N_TYPES = N_LINEAGE * N_SUB


@dataclass
class World:
    G: int
    loadings: np.ndarray        # G x Z
    type_z: np.ndarray          # N_TYPES x Z  (mean latent state per type)
    z_sd: float
    cognate: np.ndarray         # P  (gene index of each marker's partner, -1 = none)
    rho: np.ndarray             # P coupling
    prot_readout: np.ndarray    # P x Z  (protein-specific readout of latent state)
    gene_names: list
    marker_names: list


def make_world(scenario: str, G=300, P=16, Z=24, seed=0) -> World:
    rng = np.random.default_rng(seed)
    # latent programmes: 4 lineage programmes, 12 subtype programmes, 8 shared "state" programmes
    type_z = np.zeros((N_TYPES, Z))
    for t in range(N_TYPES):
        lin, sub = divmod(t, N_SUB)
        type_z[t, lin] = 3.0                         # lineage: strong
        type_z[t, 4 + t] = 1.6                       # subtype: weaker
    type_z[:, 16:] = rng.normal(0, 0.6, (N_TYPES, Z - 16))   # shared state programmes, type-dependent means
    # sparse gene loadings: each gene loads on 1-2 programmes
    loadings = np.zeros((G, Z))
    for g in range(G):
        for k in rng.choice(Z, size=rng.integers(1, 3), replace=False):
            loadings[g, k] = rng.uniform(0.4, 1.0)
    # cognate genes: markers 0-7 lineage markers (genes on lineage progs), 8-15 subtype/state markers
    cognate = np.full(P, -1)
    used = set()
    for j in range(P):
        if j < 8:
            prog = j % 4                         # markers 0-7: lineage markers (2 per lineage)
        elif j < 12:
            prog = 4 + 3 * (j - 8)               # markers 8-11: one subtype marker per lineage
        elif j < 16:
            prog = 16 + (j - 12)                 # markers 12-15: shared state programmes
        else:
            prog = int(rng.integers(0, Z))       # larger panels: extra markers on random programmes
        cand = [g for g in np.argsort(-loadings[:, prog]) if g not in used][:3]
        cognate[j] = cand[0]
        used.add(cand[0])
    prot_readout = rng.normal(0, 1, (P, Z)) * (rng.random((P, Z)) < 0.25)
    prot_readout[np.arange(P), np.argmax(np.abs(loadings[cognate]), 1)] += 0.5
    if scenario.startswith("coupled"):
        rho = np.full(P, 0.9)
    elif scenario == "decoupled":
        rho = np.full(P, 0.15)
    elif scenario == "realistic":
        rho = np.r_[np.full(8, 0.85), np.full(P - 8, 0.3)]
        cognate[[3, 7, 13, 15]] = -1                # no gene partner (DNA stain, pan-keratin, phospho, ...)
    else:
        raise ValueError(scenario)
    return World(G, loadings, type_z, 0.7, cognate, rho, prot_readout,
                 [f"g{g:03d}" for g in range(G)], [f"m{j:02d}" for j in range(P)])


def _latent(w: World, n, rng, type_probs=None):
    p = np.full(N_TYPES, 1 / N_TYPES) if type_probs is None else type_probs
    y = rng.choice(N_TYPES, size=n, p=p)
    z = w.type_z[y] + rng.normal(0, w.z_sd, (n, w.type_z.shape[1]))
    return y, z


def sample_rna(w: World, n, rng, batch_sd=0.3, depth=200):
    """Sparse RNA counts (Poisson) for n cells; returns labels, latent, log1p-normalised matrix."""
    y, z = _latent(w, n, rng)
    mu = np.exp(0.5 * (z @ w.loadings.T))
    mu *= np.exp(rng.normal(0, batch_sd, w.G))[None, :]          # dataset-level gene batch effect
    mu = mu / mu.sum(1, keepdims=True) * depth * rng.lognormal(0, 0.3, (n, 1))
    counts = rng.poisson(mu)
    x = np.log1p(counts / np.maximum(counts.sum(1, keepdims=True), 1) * 100)
    return y, z, x


def sample_protein(w: World, n, rng, n_patients=6, noise=0.5):
    """Protein intensities for n cells across patients (patient = split unit)."""
    y, z = _latent(w, n, rng)
    gene_true = 0.5 * (z @ w.loadings.T)                           # noise-free log-expression
    own = np.zeros((n, len(w.rho)))
    for j, g in enumerate(w.cognate):
        own[:, j] = gene_true[:, g] if g >= 0 else 0.0
    specific = z @ w.prot_readout.T
    own = (own - own.mean(0)) / (own.std(0) + 1e-9)
    specific = (specific - specific.mean(0)) / (specific.std(0) + 1e-9)
    rho = np.where(w.cognate >= 0, w.rho, 0.0)
    prot = rho * own + (1 - rho) * specific + rng.normal(0, noise, own.shape)
    patient = rng.integers(0, n_patients, n)
    scale = rng.lognormal(0, 0.15, (n_patients, len(rho)))
    shift = rng.normal(0, 0.2, (n_patients, len(rho)))
    prot = prot * scale[patient] + shift[patient]
    return y, z, prot, patient
