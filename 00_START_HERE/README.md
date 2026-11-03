# 00 · Start here

The repository is ordered the way the project should be read and done: **concept → literature → data → assumptions → ideas → experiment design → uncertainty → code → results.**

| Step | Folder | Read it to… | Time |
|---|---|---|---|
| 0 | [`00_START_HERE/project_concept.md`](project_concept.md) | understand the question, Justin's H1 vs H2, and the information-map layer | 10 min |
| 1 | [`01_literature_review/`](../01_literature_review/) | know each algorithm's idea, **input format**, output and role (6 families) | 1–2 h |
| 2 | [`02_datasets/`](../02_datasets/) | pick data by **RNA–protein pairing tier**; [`dataset_playbook.md`](../02_datasets/dataset_playbook.md) says what to do with each | 30 min |
| 3 | [`03_assumptions/`](../03_assumptions/) | the assumptions behind any "informative/transferable" claim, and **why information fails to transfer** | 30 min |
| 4 | [`04_ideas/`](../04_ideas/) | one folder per idea (I1, H1, H2, N2–N7) with status, data, methods, UQ, cost | 20 min |
| 5 | [`05_experiment_design/`](../05_experiment_design/) | E1–E8, evaluation protocol, transferability diagnostics | 40 min |
| 6 | [`06_uncertainty_quantification/`](../06_uncertainty_quantification/) | bootstrap, permutation, conformal, model uncertainty, multiplicity | 30 min |
| 7 | [`07_code/`](../07_code/) | run the information map + UQ (runnable now); H1 protein foundation-model template; H2 stub | – |
| 8 | [`08_results/`](../08_results/) | simulation self-checks; real-data results go here | – |

Glossary: [`glossary.md`](glossary.md).

## The 4-week plan (suggested)

| Week | Do | Output |
|---|---|---|
| 1 | Xenium Protein RCC (P1) → E1–E4 + UQ; H1 go/no-go (both models embed SPATCH colon) | first real information map |
| 2 | replicate on P2/P6; E5 region hold-out; E6 residuals; H2 stage 1 | transfer panel, residual spatial structure |
| 3 | E7 (P1 → DBiTplus / spatial-CITE-seq); E8 + SPATCH registration; H2 stage 2 or H1 fine-tuning | transfer and pairing results |
| 4 | test the pre-registered prediction (H1/H2 gains where the map says); figures; report | final report |
