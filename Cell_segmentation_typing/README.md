### Cell segmentation
- The routine segmentation pipeline is based [ClusterMap package](https://github.com/wanglab-broad/ClusterMap)
- Example scripts for cell segmentation in sagittal sections are provided

### Combined cell typing
- Cells in all samples are pooled together for combined analysis
- The analysis includes the following steps:
  - Preprocessing
  - Batch effect check and removal
  - Integration with [scRNA-seq data](http://mousebrain.org/adolescent/downloads.html)  
  - Leiden Clustering
  - Label transfer of main clusters

### Subcluster cell types
- We further subcluster by Leiden clustering and then manually annotate within each main cluster
- Example script of 'Di- and mesencephalon excitatory neurons' and 'cerebellum neurons' is provided


### Plot UMAP and cell type map
- [UMAP](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.umap.html#scanpy-tl-umap) is plotted colored by subclusters.
- Spatial cell-type maps are plotted colored by subclusters of 20 mouse CNS slices. Each dot represents one cell.
