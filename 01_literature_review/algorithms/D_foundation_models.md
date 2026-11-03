# D. Foundation models (H1 transfer, N7 probes)

**Use when** you want knowledge from a large pretraining corpus of one modality to help a small task in the other modality.
**Key question for us:** what *carries* the transfer? Shared feature identity (gene ↔ protein token), shared tissue structure (graph), or nothing?

---

## D1. Nicheformer (Tejada-Lapuerta, Schaar et al., *Nat Methods* 2025): RNA side

**Pretraining.** SpatialCorpus-110M: about 57M dissociated and 53M spatially resolved cells, human and mouse. Vocabulary of 20,310 gene tokens (orthologous protein-coding genes plus species-specific genes). Transformer with a 1,500-token context.

**Input format (from `theislab/nicheformer`).**
1. AnnData with **raw counts**, genes reordered to the model's vocabulary (template: `data/model_means/model.h5ad`).
2. Library-size normalisation (`sf_normalize`), then division by the **technology mean** vector (`xenium_mean_script.npy`, `cosmx_…`, `merfish_…`, `iss_…`, `dissociated_…`). This puts different assays on a comparable scale.
3. Each cell becomes the **ranked list of non-zero genes** (highest relative expression first) → integer token IDs (offset by 30 reserved tokens), padded or truncated.
4. Contextual tokens: species, assay, modality.

Tokenisation notebooks: `notebooks/tokenization/xenium_human_colon.ipynb`, `cosmx_human_lung.ipynb`, …; pretrained weights on Mendeley Data.

**Feeding protein into Nicheformer (H1, protein → model).** Map each marker to its HGNC gene (MaxFuse table), drop unmapped markers (or add new tokens and fine-tune), rank the cell's markers by arcsinh-scaled intensity, and tokenise. Only ~16–60 tokens per cell are non-zero, far fewer than for RNA. That **panel-overlap** effect is ablation A3 in Justin's plan.

---

## D2. Spatium (bioRxiv 2026, doi 10.64898/2026.07.23.740264): protein side

**Pretraining (repo README).** About 51M cells, 70 datasets, 15 tissues, 18 diseases, 4 technologies (IMC, MIBI, CODEX, CyCIF), protein vocabulary of 234. A transformer over protein tokens ordered by within-cell relative abundance, with metadata tokens:
```
[CLS] [Technology] [Tissue] [Health/Disease] [Disease Type] [Protein 1] [Protein 2] ... [PAD]
```
Self-supervised by recovering masked proteins. It supports zero-shot reconstruction of an unobserved protein from the others.

**Input format (from `ploughhh/Spatium`).**
```python
schema = DataSchema(); schema.set(adata, Technology="CODEX", Tissue_type="Intestine",
                                  Health_status="Health", Disease_type="Normal")
adata.write_h5ad(path)                       # tokenizer reads from disk
tok = FineTuneTokenizer(all_proteins=all_proteins, protein_id_map=protein_id_map,
                        side_maps=side_maps, protein_token_base=PROTEIN_TOKEN_BASE)
tok.collect_data([path], cell_type_col="Cell Type"); tokens = tok.tokenize()   # -> zarr: tokens [N, L], label [N]
```
`src/protein_standard_mapping.tsv` and `protein_gene_conversion.csv` standardise antibody names. A separate **spatial niche tokenizer** handles neighbourhood tasks.

**Feeding RNA into Spatium (H1, RNA → model).** Keep only genes that code for one of the 234 vocabulary proteins, rank by normalised expression, and tokenise as protein tokens. The input therefore shrinks to the mappable genes.

---

## D3. HEIST (Madhu et al., ICLR 2026; arXiv 2506.11152)

**Idea.** A *hierarchical graph* transformer. The upper level is a spatial cell graph built from **Voronoi adjacency**. Each cell links to a lower-level **gene co-expression network for its cell type**: genes are MAGIC-denoised and pairs whose mutual information exceeds τ are connected. Cross-level message passing lets the microenvironment reshape gene embeddings and vice versa. Gene embeddings are computed from the network rather than a fixed vocabulary, so it can be applied to **unseen features, including proteins**.

**Pretraining.** 22.3M cells, 124 tissues, 15 organs; spatially aware contrastive plus masked-autoencoding objectives.

**Input.** Normalised HVG matrix, coordinates, cell-type labels (or Leiden clusters). For protein data the "genes" are markers, and co-expression networks are built per cell type from marker MI.

**Relevance.** The only published model **pretrained on SPT and evaluated on protein (CODEX, MIBI)**. Per Justin's review, frozen embeddings + MLP beat the best baseline on 5 of 6 outcome tasks but lacked from-scratch and protein-pretrained comparators. Our H1 contrasts C1/C2 supply exactly those missing controls.

---

## D4. Novae (*Nat Methods* 2025; `MICS-Lab/novae`)

A graph-based foundation model for **spatial domains** in SPT, working across panels and technologies (zero-shot domain assignment):
```python
novae.spatial_neighbors(adata)
model = novae.Novae.from_pretrained("prism-oncology/novae-human-0")
model.compute_representations(adata, zero_shot=True); model.assign_domains(adata)
```
**Use.** H1 fallback if Nicheformer will not run. It also provides region labels for E5 (region hold-out) that are defined **from RNA only**, which avoids circularity when testing protein transfer.

## D5. KRONOS (Mahmood lab; arXiv 2506.03373) and VirTues (Bunne lab; *Nature* 2026)

Image-based protein foundation models.
- **KRONOS:** self-supervised on 47M single-marker patches covering 175 markers, 16 tissues and 8 platforms; segmentation-free patch embeddings; weights on Hugging Face (gated, academic use). KRONOS2 has been announced.
- **VirTues:** marker-aware multi-scale transformer; marker tokens come from **protein language model** embeddings (one `[UniprotID].pt` per marker), so it handles arbitrary panels.

**Use.** Protein-side baselines for H1/H2 (the student must be competitive with them). Neither sees RNA, so they bound what protein alone achieves.

---

## How to compare foundation models fairly (from H1, generalised)

| Control | Why |
|---|---|
| Same architecture, random init | separates "pretraining helps" from "this architecture is good" |
| In-modality pretrained model | ceiling for transfer efficiency |
| PCA + logistic regression / Leiden | foundation models often barely beat it (Kedzierska et al., *Genome Biol* 2025: zero-shot scGPT/Geneformer trailed HVG and scVI) |
| Pretraining-leakage check | the test slides may be in the pretraining corpus |
| Token-identity ablations | A1 re-init mapped tokens; A2 random token mapping; A3 shrink panel overlap |
