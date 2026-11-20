# A2. Feature identity: a protein marker and its gene measure the same molecular entity

**Statement.** The antibody for marker *m* binds the product of gene *g(m)*, and the cognate map is correct.

**Typical violations.**
- **Isoforms:** CD45RA and CD45RO are both *PTPRC* (Xenium added isoform-specific probes for exactly this reason).
- **Multi-gene antibodies:** pan-cytokeratin (many *KRT* genes), HLA-DR (*HLA-DRA/DRB1/DRB5*).
- **Post-translational states:** phospho-S6, cleaved caspase-3, Ki-67 vs *MKI67* timing. The gene is right but the state is not encoded in mRNA.
- **No gene at all:** DNA stains, collagen autofluorescence, some carbohydrate epitopes.
- **Clone differences** across cohorts: the same marker name can be a different measurement.

**If violated.** E2 underestimates cognate information, and H1 token mapping (gene token ↔ protein token) transfers the wrong prior.

**How we test it.**
- Keep a **mapping-confidence label** for each marker: 1:1, 1:many, isoform, PTM, none. Report E1/E2 results stratified by that label.
- H1 ablation A2 (map markers to random gene tokens) and A1 (re-initialise mapped tokens).
- For 1:many markers, compare "best single gene" vs "sum of family genes" as the cognate input.
