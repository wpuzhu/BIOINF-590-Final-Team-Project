"""
Torch / PyTorch Geometric models for Stage 2 (cluster only; not importable without torch).

Student : protein of the cell + 2-layer GraphSAGE over the radius graph  -> fine cell type
Teacher : the student's encoder + an RNA encoder (cell RNA PCs, or neighbourhood RNA for
          adjacent sections), late fusion                                -> fine cell type
Heads   : classifier; auxiliary projection head (feature KD / RNA reconstruction / CCA target);
          contrastive projection head (InfoNCE between student embedding and RNA embedding)
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv


class GraphEncoder(nn.Module):
    def __init__(self, d_in, hidden=128, use_graph=True, dropout=0.2):
        super().__init__()
        self.use_graph = use_graph
        self.inp = nn.Linear(d_in, hidden)
        if use_graph:
            self.c1, self.c2 = SAGEConv(hidden, hidden), SAGEConv(hidden, hidden)
        else:
            self.c1, self.c2 = nn.Linear(hidden, hidden), nn.Linear(hidden, hidden)
        self.drop = dropout
        self.out_dim = hidden

    def forward(self, x, edge_index):
        h = F.relu(self.inp(x))
        if self.use_graph:
            h = F.relu(self.c1(h, edge_index)); h = F.dropout(h, self.drop, self.training)
            h = F.relu(self.c2(h, edge_index))
        else:
            h = F.relu(self.c1(h)); h = F.dropout(h, self.drop, self.training)
            h = F.relu(self.c2(h))
        return h


class Student(nn.Module):
    def __init__(self, d_prot, n_classes, hidden=128, use_graph=True, aux_dim=0, proj_dim=64):
        super().__init__()
        self.enc = GraphEncoder(d_prot, hidden, use_graph)
        self.cls = nn.Linear(hidden, n_classes)
        self.aux = nn.Linear(hidden, aux_dim) if aux_dim else None
        self.proj = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, proj_dim))

    def forward(self, x, edge_index):
        h = self.enc(x, edge_index)
        return dict(h=h, logits=self.cls(h), aux=self.aux(h) if self.aux is not None else None, z=self.proj(h))


class Teacher(nn.Module):
    def __init__(self, d_prot, d_rna, n_classes, hidden=128, use_graph=True):
        super().__init__()
        self.enc = GraphEncoder(d_prot, hidden, use_graph)
        self.rna = nn.Sequential(nn.Linear(d_rna, hidden), nn.ReLU(), nn.Dropout(0.2), nn.Linear(hidden, hidden), nn.ReLU())
        self.fuse = nn.Linear(2 * hidden, hidden)
        self.cls = nn.Linear(hidden, n_classes)

    def forward(self, x, edge_index, r):
        h = F.relu(self.fuse(torch.cat([self.enc(x, edge_index), self.rna(r)], 1)))
        return dict(h=h, logits=self.cls(h))


class RNAProjector(nn.Module):
    """RNA side of the contrastive alignment (trained jointly with the student)."""
    def __init__(self, d_rna, hidden=128, proj_dim=64):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_rna, hidden), nn.ReLU(), nn.Linear(hidden, proj_dim))

    def forward(self, r):
        return self.net(r)


def info_nce(z_a, z_b, tau=0.1, max_n=4096):
    """Symmetric InfoNCE on a random subset (memory): positives are the same cell."""
    n = z_a.shape[0]
    if n > max_n:
        idx = torch.randperm(n, device=z_a.device)[:max_n]
        z_a, z_b = z_a[idx], z_b[idx]
    z_a, z_b = F.normalize(z_a, dim=1), F.normalize(z_b, dim=1)
    logits = z_a @ z_b.T / tau
    y = torch.arange(len(z_a), device=z_a.device)
    return 0.5 * (F.cross_entropy(logits, y) + F.cross_entropy(logits.T, y))


def kd_loss(logits_s, p_teacher_T, T):
    return F.kl_div(F.log_softmax(logits_s / T, 1), p_teacher_T, reduction="batchmean") * T * T
