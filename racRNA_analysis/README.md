### RNA barcode analysis
Inputs:
- racRNA spot data and mRNA spot data (with cell labels)
- `pd_tissue.csv` file containing Rank4/5/Tissue annotations for each cell
- `color_final.csv` file specifying order/color of cell types and regions in the plots

Files:
- `racRNA_classification.ipynb` assigns racRNA spots into cells, outputs two files for each sample:
-- outputs a file of each racRNA spot and its associated cell
-- outputs a file containing the number of racRNA in each cell with nonzero racRNA count
-- currently configured to run on Well 11

- `figure_generation.ipynb` creates figure using the number of racRNA in each cell as input, producing 6 files:
-- outputs two files containing sagittal and coronal input data for figures
-- creates two boxplots for coronal (linear scaled x-axis and logarithmically scaled x-axis) and two boxplots for sagittal
-- currently configured to produce figures based on Rank 5 cell type (Figure 4d)
