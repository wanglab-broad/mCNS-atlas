### Intermediate mapping
- Specifically, for each of all 1022 genes in the STARmap PLUS, we performed an intermediate mapping to align each STARmap PLUS cell with the most similar set of cells in the scRNA-seq dataset: `intermediate_mapping.py`

### Final imputation
- We first checked the quality of the scRNA-seq data1: genes with average read < 0.005 / sum read < 740 across 146,201 cells (50th percentile of the data) were filtered; genes with maximum read <= 10 were filtered. 
- Then, for every cell, we calculated each gene’s imputed expression level as the weighted average of the gene’s expression across the associated set of scRNA-seq atlas cells, where weights were proportional to the number of times each scRNA-seq atlas cell was present: `final_imputation_compute.py`
