"""
A small numpy MLP that supports every Stage 2 mechanism in one loss, so that the
simulation compares mechanisms with identical architecture and optimiser:

    L = lam * CE(y_labelled)                               hard labels (labelled cells only)
      + (1 - lam) * T^2 * KL(p_teacher^T || p_student^T)   distillation (all training cells)
      + beta * ||head(h) - target||^2                      auxiliary target (all training cells)

`target` is what makes the mechanism: the teacher's hidden layer (feature distillation), the
RNA principal components (auxiliary RNA reconstruction) or the RNA part shared with protein
(CCA variates; a closed-form stand-in for contrastive alignment). The torch pipeline uses the
same loss with InfoNCE for alignment (stage2_privileged_transfer/pipeline).
"""
from __future__ import annotations

import numpy as np


def softmax(z):
    z = z - z.max(1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(1, keepdims=True)


class NumpyMLP:
    def __init__(self, n_classes, hidden=64, lr=1e-2, epochs=300, wd=1e-4, seed=0):
        self.K, self.h, self.lr, self.epochs, self.wd, self.seed = n_classes, hidden, lr, epochs, wd, seed

    def _init(self, d, m):
        r = np.random.default_rng(self.seed)
        self.W1 = r.normal(0, np.sqrt(2 / d), (d, self.h)); self.b1 = np.zeros(self.h)
        self.W2 = r.normal(0, np.sqrt(1 / self.h), (self.h, self.K)); self.b2 = np.zeros(self.K)
        self.Wa = r.normal(0, np.sqrt(1 / self.h), (self.h, max(m, 1))); self.ba = np.zeros(max(m, 1))

    def fit(self, X, y=None, labelled=None, soft=None, T=2.0, lam=1.0, target=None, beta=0.0,
            smooth=0.0):
        """X all training cells; y labels with `labelled` boolean mask; soft teacher probs at
        temperature T for all cells (or None); target aux matrix (or None)."""
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-6
        X = (X - self.mu) / self.sd
        n, d = X.shape
        m = 0 if target is None else target.shape[1]
        if target is not None:
            tmu, tsd = target.mean(0), target.std(0) + 1e-6
            target = (target - tmu) / tsd
        self._init(d, m)
        params = ["W1", "b1", "W2", "b2", "Wa", "ba"]
        mom = {k: np.zeros_like(getattr(self, k)) for k in params}
        vel = {k: np.zeros_like(getattr(self, k)) for k in params}
        L = np.where(labelled)[0] if labelled is not None else np.arange(0)
        Y = None
        if len(L):
            Y = np.full((len(L), self.K), smooth / self.K)
            Y[np.arange(len(L)), y[L]] += 1 - smooth
        for t in range(1, self.epochs + 1):
            Z1 = X @ self.W1 + self.b1
            H = np.maximum(Z1, 0)
            Z = H @ self.W2 + self.b2
            dZ = np.zeros_like(Z)
            if len(L) and lam > 0:
                dZ[L] += lam * (softmax(Z[L]) - Y) / len(L)
            if soft is not None and lam < 1:
                dZ += (1 - lam) * T * (softmax(Z / T) - soft) / n
            dH = dZ @ self.W2.T
            g = {}
            if target is not None and beta > 0:
                A = H @ self.Wa + self.ba
                dA = 2 * beta * (A - target) / (n * m)
                dH += dA @ self.Wa.T
                g["Wa"], g["ba"] = H.T @ dA, dA.sum(0)
            else:
                g["Wa"], g["ba"] = np.zeros_like(self.Wa), np.zeros_like(self.ba)
            dZ1 = dH * (Z1 > 0)
            g["W2"], g["b2"] = H.T @ dZ, dZ.sum(0)
            g["W1"], g["b1"] = X.T @ dZ1, dZ1.sum(0)
            for k in params:
                gk = g[k] + (self.wd * getattr(self, k) if k.startswith("W") else 0)
                mom[k] = 0.9 * mom[k] + 0.1 * gk
                vel[k] = 0.999 * vel[k] + 0.001 * gk ** 2
                mh, vh = mom[k] / (1 - 0.9 ** t), vel[k] / (1 - 0.999 ** t)
                setattr(self, k, getattr(self, k) - self.lr * mh / (np.sqrt(vh) + 1e-8))
        return self

    def hidden(self, X):
        return np.maximum(((X - self.mu) / self.sd) @ self.W1 + self.b1, 0)

    def predict_proba(self, X, T=1.0):
        return softmax((self.hidden(X) @ self.W2 + self.b2) / T)
