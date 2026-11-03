# 04 · Ideas

One folder per idea. Each README has the same header: status, question, assumptions, data, methods, experiments, UQ and cost. Ideas are listed from **foundational and cheap** to **ambitious**.

| Order | Idea | One line | Depends on |
|---|---|---|---|
| 1 | [**I1 information map**](I1_information_map/) | which RNA inputs inform each protein, at what scale, and whether that transfers | – (runnable now) |
| 2 | [N2 residual ctSVPs](N2_residual_ctSVP/) | protein spatial structure left after RNA is accounted for | I1 predictor |
| 3 | [N1 spatial-variability transfer](N1_spatial_variability_transfer/) | do RNA ctSVGs become protein ctSVPs? | v1 detectors |
| 4 | [N3 minimal gene panel](N3_minimal_gene_panel/) | which transcripts suffice to predict the protein panel | I1 |
| 5 | [**H2 privileged information**](H2_privileged_information/) | RNA teacher → protein-only student | I1 (where to expect gains), E8 |
| 6 | [**H1 transfer learning**](H1_transfer_learning/) | Nicheformer ↔ Spatium cross-modal pretraining | I1, N7 |
| 7 | [N7 token-identity probe](N7_token_identity_probe/) | do pretrained embeddings encode the other modality? | H1 embeddings |
| 8 | [N4 virtual protein](N4_virtual_protein/) | trusted protein predictions with intervals for RNA-only cohorts | I1, UQ |
| 9 | [N5 shared vs specific factors](N5_shared_vs_specific_factors/) | joint factor decomposition of RNA and protein | – |
| 10 | [N6 counterfactual neighbourhoods](N6_counterfactual_neighbourhoods/) | what happens to protein if context changes | I1 E3/E4 |

**Recommended scope for 4 weeks:** I1 + N2 as the guaranteed core, plus **one** of H2 (most novel; code exists in simulation) or H1 (needs the Week-1 go/no-go). That matches option C in Justin's "H1 vs H2" deck (H2 plus one H1 arm) with I1 underneath.
