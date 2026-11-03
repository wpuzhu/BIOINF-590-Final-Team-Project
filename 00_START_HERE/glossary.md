# Glossary

| Term | Meaning here |
|---|---|
| **SPT** | spatial transcriptomics (Xenium, CosMx, MERFISH, Visium/Visium HD, Stereo-seq, DBiT-seq) |
| **Spatial proteomics** | multiplexed protein imaging (CODEX/PhenoCycler, IMC, MIBI, CyCIF) or sequencing-based antibody tags (ADT: spatial-CITE-seq, SPOTS, Visium protein) |
| **Cognate gene** | the gene encoding a protein marker (CD3 → *CD3E*). The mapping can be 1:1, 1:many, isoform, PTM (phospho) or none |
| **Unit** | the independent thing we split, resample and permute: patient, slide, or spatial block. Never a cell |
| **Information map** | per protein × cell type: held-out variance explained by cell type, cognate gene, other genes and neighbourhood RNA (idea I1) |
| **Informative** | adding an input raises held-out R² in units from the *same* population (CI excludes 0) |
| **Transferable** | the gain persists in units from a *different* population: region, dataset, platform, or imperfect pairing |
| **Composition vs arrangement** | composition: which cell types are nearby; arrangement: what state the nearby cells are in. Separated by E4 shuffles |
| **Covariate shift** | target has a different mix of inputs (cell types, contexts); P(protein \| RNA) unchanged |
| **Concept shift** | P(protein \| RNA) itself differs in the target |
| **Support shift** | target contains inputs never seen in training (extrapolation) |
| **Residual ctSVP** | protein spatial variability within a cell type that remains after removing the RNA-predicted part |
| **Cluster bootstrap** | resample units (not cells) with replacement |
| **Sign-flip test** | paired permutation test over units: flip the sign of each unit's difference |
| **Toroidal shift** | spatial null that translates one field with wrap-around, keeping its autocorrelation |
| **Mondrian conformal** | conformal intervals calibrated separately per group (here: per cell type) |
| **Weighted conformal** | conformal with calibration residuals reweighted by the density ratio target/source (covariate shift) |
| **Privileged information (H2)** | data available only during training (RNA), never at test time |
| **Teacher / student (H2)** | teacher sees protein + RNA; student sees protein only and learns from labels + teacher outputs |
| **C1–C4 (H1)** | pre-registered cross-modal contrasts in Justin's H1 plan |
| **A1–A8** | assumptions in `03_assumptions/` |
| **E1–E8** | experiments in `05_experiment_design/` |
| **U1–U5** | UQ method families in `06_uncertainty_quantification/` |
