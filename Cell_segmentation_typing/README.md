### cell segmentation
- The segmentation is based [ClusterMap package](https://github.com/wanglab-broad/ClusterMap)
- Example scripts for cell segmentation in sagittal sections are provided

### combined cell typing
- Cells in all samples are pooled together for combined analysis
- The analysis includes the following steps:
  - Preprocessing
  - Batch effect check and removal
  - Integration with [scRNA-seq data](http://mousebrain.org/adolescent/downloads.html)  
  - Leiden Clustering
  - Label transfer
