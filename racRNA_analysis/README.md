### RNA barcode analysis
- Assign circular RNA barcode spots into cells

After spot-calling of circular RNA barcode spots, we stitched the spots in each tile together based on tile location information and then assigned circular RNA barcode spots into cells identified by endogenous genes. In order to include only those spots with high confidence of being in a cell, we used a Nearest Neighbors algorithm to determine which RNA barcode amplicons are in which cells. We first used sklearn.neighbors.NearestNeighbors to identify the mRNA spots closest to each RNA barcode spot. We considered a circular RNA barcode spot to be in a cell if its closest mRNA spot was within a distance of 3.5 pixels, or if its three closest mRNA spots were all within a distance of 10 pixels and were all in the same cell. Circular RNAs not satisfying one of the two conditions were considered to not be in any cells. We then counted the total number of circular RNA barcodes for each cell.

- Cell type-based statistics

We first grouped the cells by subtypes. Then, for each cell type, we computed summary statistics of the 2.5th, 25th, 50th, 75th, and 97.5th percentiles using pandas.core.groupby.DataFrameGroupBy.quantile to generate a boxplot of circular RNA barcode expression by cell type in both coronal and sagittal samples.

- Tissue region-based statistics

We similarly computed the 2.5th, 25th, 50th, 75th, and 97.5th percentiles for each tissue region after grouping cells by the tissue regions as generated above. 
