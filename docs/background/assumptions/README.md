# 03 · Assumptions

Every claim of the form "RNA information X *informs* / *transfers to* protein Y" rests on the assumptions below. Each file states the assumption, how it fails, what the failure would look like in our results, and which experiment or UQ check tests it.

| ID | Assumption | Tested by |
|---|---|---|
| [A1](A1_cell_correspondence.md) | RNA and protein in a unit come from the same cells | E8 jitter/pairing sweep; MaxFuse vs registration |
| [A2](A2_feature_identity.md) | a marker and its cognate gene measure the same entity | mapping-confidence strata; H1 A1/A2 ablations |
| [A3](A3_label_comparability.md) | cell-type labels are comparable and not circular | RNA-derived vs protein-derived labels |
| [A4](A4_stationarity.md) | P(protein \| RNA) is the same in source and target | E5, E7, shift decomposition, weighted conformal |
| [A5](A5_exchangeable_units.md) | resampling units are independent | block size vs autocorrelation range; toroidal vs naive nulls |
| [A6](A6_spatial_scale.md) | the neighbourhood radius matches the biology | E3 sweep with bootstrap CI on the peak |
| [A7](A7_measurement_models.md) | noise models / transforms suit each modality | isotype controls; spillover check |
| [A8](A8_no_leakage.md) | nothing about test units leaks into training | pretraining corpus check; nested CV |

**Why information fails to transfer** (biological vs technical mechanisms, their signatures, and a decision tree): [`why_information_does_not_transfer.md`](why_information_does_not_transfer.md).
