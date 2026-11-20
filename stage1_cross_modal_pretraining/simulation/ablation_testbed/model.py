"""
A small "foundation model" with the two kinds of pretrained knowledge the H1 ablations separate:

  layer 1 = TOKEN EMBEDDING TABLE  (vocabulary x d): what each gene token means   -> "feature identity"
  layer 2 = SHARED ENCODER          (d x d):          how token evidence combines  -> "structure / capacity"

Cells are tokenised like Nicheformer/Spatium: features ranked within the cell, the top-L get
weights that decay with rank, placed at their token's column in the vocabulary. The network
is an embedding-bag MLP (sklearn), pretrained by masked reconstruction on an RNA corpus and
then fine-tuned (or frozen + probed) on a protein cell-typing task.

Arms (what is handed to the protein task):
  scratch        random init everywhere                                     (H1 control C1)
  noise_pre      pretrained on an RNA corpus with genes shuffled independently across cells:
                 same marginals, no co-expression -> generic optimisation benefit only
  pretrained     RNA-pretrained, protein markers mapped to their cognate gene tokens
  A1_reinit      pretrained encoder kept, embedding rows of the mapped tokens re-initialised
  A2_random_map  pretrained, each marker mapped to a RANDOM (wrong) gene token
  A3_overlap_q   pretrained, only a fraction q of markers mapped to cognate tokens; the rest get
                 NEW tokens with fresh embeddings (information content unchanged)
"""
from __future__ import annotations

import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPClassifier, MLPRegressor

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", message="Got `batch_size` less than 1")   # 5 labels x 12 types < batch 64: sklearn clips, harmless
HIDDEN = (64, 64)


def tokenize(X, cols, V, L=40):
    """Rank features within each cell; top-L get weight 1 - r/L at their vocabulary column."""
    n = X.shape[0]
    T = np.zeros((n, V))
    order = np.argsort(-X, axis=1)
    for r in range(min(L, X.shape[1])):
        f = order[:, r]
        val = X[np.arange(n), f]
        w = (1 - r / L) * (val > -np.inf)
        T[np.arange(n), np.asarray(cols)[f]] = np.where(np.isfinite(val), w, 0)
    return T


def rna_tokens(x, gene_mean, V, L=40):
    rel = x / (gene_mean + 1e-6)
    rel[x <= 0] = -np.inf                               # only expressed genes are tokens
    return tokenize(rel, np.arange(x.shape[1]), V, L)


def protein_tokens(prot, patient, cols, V):
    z = prot.copy()
    for p in np.unique(patient):                        # per-patient z-score, then rank within cell
        m = patient == p
        z[m] = (z[m] - z[m].mean(0)) / (z[m].std(0) + 1e-9)
    return tokenize(z, cols, V, L=prot.shape[1])


# ------------------------------------------------------------------ pretraining
def pretrain(T_corpus, x_corpus, seed=0, mask_frac=0.3, n_comp=30, max_iter=80):
    """Self-supervised pretraining on the RNA corpus (no labels).

    Input : the cell's rank tokens with a random 30% hidden (masked-token style).
    Target: the cell's DENOISED full expression profile, i.e. its projection on the top-30
            principal components of the corpus (unsupervised low-rank structure, MAGIC/PCA-like).
    To hit the target from a partial, noisy token set, the network has to learn which genes
    co-vary (embeddings) and how programmes combine into cell states (encoder): the two kinds
    of knowledge the ablations try to separate. (A raw masked-token target was tried first; it
    is dominated by Poisson noise and learned nothing transferable - see README.)"""
    xs = (x_corpus - x_corpus.mean(0)) / (x_corpus.std(0) + 1e-6)
    pca = PCA(n_comp, random_state=seed).fit(xs)
    Y = pca.inverse_transform(pca.transform(xs))
    rng = np.random.default_rng(seed)
    M = rng.random(T_corpus.shape) < mask_frac
    reg = MLPRegressor(hidden_layer_sizes=HIDDEN, activation="relu", alpha=1e-4, batch_size=256,
                       learning_rate_init=1e-3, max_iter=max_iter, tol=1e-6, n_iter_no_change=30,
                       random_state=seed)
    reg.fit(T_corpus * (~M), Y)
    return [c.copy() for c in reg.coefs_[:2]], [b.copy() for b in reg.intercepts_[:2]], reg.loss_curve_


def shuffle_genes(x, seed=0):
    """Permute each gene independently across cells: marginals kept, co-expression destroyed."""
    rng = np.random.default_rng(seed)
    xs = x.copy()
    for g in range(x.shape[1]):
        xs[:, g] = xs[rng.permutation(x.shape[0]), g]
    return xs


# ------------------------------------------------------------------ arms
def column_map(arm, cognate, G, rng, q=1.0):
    """Vocabulary column for each marker under an arm, plus the set of columns to re-initialise."""
    P = len(cognate)
    cols = np.where(cognate >= 0, cognate, G + np.arange(P))          # no partner -> new token
    reinit = set(G + np.arange(P))                                     # new tokens always fresh
    mappable = np.where(cognate >= 0)[0]
    if arm == "A2_random_map":
        pool = np.setdiff1d(np.arange(G), cognate[mappable])
        cols[mappable] = rng.choice(pool, size=len(mappable), replace=False)
    elif arm.startswith("A3_overlap"):
        keep = rng.choice(mappable, size=int(round(q * len(mappable))), replace=False)
        drop = np.setdiff1d(mappable, keep)
        cols[drop] = G + drop
    if arm == "A1_reinit":
        reinit |= set(cols.tolist())
    return cols, sorted(reinit)


def _glorot(rng, fan_in, fan_out, shape):
    lim = np.sqrt(6 / (fan_in + fan_out))
    return rng.uniform(-lim, lim, shape)


def init_weights(pre, reinit_cols, V, rng):
    if pre is None:
        return None
    (W1, W2), (b1, b2) = pre
    W1 = W1.copy()
    for c in reinit_cols:
        W1[c] = _glorot(rng, V, HIDDEN[0], HIDDEN[0])
    return [W1, W2.copy()], [b1.copy(), b2.copy()]


def finetune(Xtr, ytr, weights, n_classes, seed, max_iter=150):
    clf = MLPClassifier(hidden_layer_sizes=HIDDEN, activation="relu", alpha=1e-3, batch_size=64,
                        learning_rate_init=1e-3, max_iter=max_iter, random_state=seed)
    if weights is not None:
        clf.partial_fit(Xtr, ytr, classes=np.arange(n_classes))           # allocate layers (1 epoch)
        clf.set_params(warm_start=True)
        rng = np.random.default_rng(seed)
        clf.coefs_[0], clf.coefs_[1] = [w.copy() for w in weights[0]]
        clf.intercepts_[0], clf.intercepts_[1] = [b.copy() for b in weights[1]]
        clf.coefs_[2] = _glorot(rng, HIDDEN[1], n_classes, clf.coefs_[2].shape)
        clf.intercepts_[2] = np.zeros_like(clf.intercepts_[2])
        del clf._optimizer                                                   # fresh Adam state
    clf.fit(Xtr, ytr)
    return clf


def embed(X, weights):
    (W1, W2), (b1, b2) = weights
    h = np.maximum(X @ W1 + b1, 0)
    return np.maximum(h @ W2 + b2, 0)


def frozen_probe(Xtr, ytr, Xte, weights):
    Htr, Hte = embed(Xtr, weights), embed(Xte, weights)
    mu, sd = Htr.mean(0), Htr.std(0) + 1e-6
    lr = LogisticRegression(C=1.0, max_iter=3000)
    lr.fit((Htr - mu) / sd, ytr)
    return lr.predict((Hte - mu) / sd)


def random_weights(V, rng):
    return ([_glorot(rng, V, HIDDEN[0], (V, HIDDEN[0])), _glorot(rng, HIDDEN[0], HIDDEN[1], HIDDEN)],
            [np.zeros(HIDDEN[0]), np.zeros(HIDDEN[1])])
