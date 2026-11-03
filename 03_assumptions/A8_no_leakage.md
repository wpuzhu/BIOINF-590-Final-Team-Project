# A8. No leakage between training and evaluation

**Statement.** Nothing about the test units influences training: not preprocessing parameters, not pretraining corpora, not hyper-parameter selection, not labels.

**Common leaks and guards.**

| Leak | Guard |
|---|---|
| Foundation model pretrained on the test slides | check Nicheformer's SpatialCorpus list and Spatium's study list (H1 leakage check); report seen vs unseen separately |
| Hyper-parameters tuned on test patients | nested CV or a separate validation unit (H2 plan: validation patients choose λ, T, β) |
| Scalers/PCA fit on all cells | unsupervised, so a minor leak; fit on training units for final results (the demo fits PCA once on intact RNA for speed) |
| Neighbourhood features crossing split boundaries | with block splits, neighbour features for a test cell may use training cells' RNA. That is allowed for *RNA* (an input), but never for *protein* (the target) |
| Cell types called from protein, then used as inputs | see A3; interpret gains only, or use RNA-derived types |
