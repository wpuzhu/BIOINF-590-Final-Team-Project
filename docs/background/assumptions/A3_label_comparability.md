# A3. Cell-type labels are comparable and not circular

**Statement.** Cell types are defined the same way across samples and modalities, and the labels used to *evaluate* a protein model are not derived from the RNA that *trains* it.

**Violations.**
- Labels clustered from RNA, then used as targets for an RNA-taught protein student (H2): the teacher wins trivially.
- Different annotation granularity across datasets (e.g. "T cell" vs CD4/CD8 subsets).
- Labels from protein, then used as an **input** to predict that protein: cell-type features "leak" protein information.

**Consequences for our design.**
- In E1, the cell-type step uses labels; if they were called from protein, step 1 already contains protein information. That is fine **as long as we only interpret *gains* over step 1**. Better: also run E1 with **RNA-derived** cell types, where the cell-type step is purely RNA information.
- H2's main claim uses protein + morphology labels (Justin's plan); RNA-derived labels are reported as a separate cross-modal task.
