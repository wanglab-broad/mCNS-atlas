### Cell segmentation
- The routine segmentation pipeline is based [ClusterMap package](https://github.com/wanglab-broad/ClusterMap)
- Example scripts for cell segmentation in sagittal sections are under `segmentation`. Change input for cell segmentation of other sections.

### Combined cell typing
- Cells in all samples are pooled together for combined analysis in `combined_cell_typing.ipynb`
- The analysis includes the following steps:
  - Preprocessing
  - Batch effect check and removal
  - Integration with [scRNA-seq data](http://mousebrain.org/adolescent/downloads.html)  
  - Leiden Clustering
  - Label transfer of main clusters

### Subcluster cell types
- We further subcluster by Leiden clustering and then manually annotate within each main cluster
- Example script of 'Di- and mesencephalon excitatory neurons' and 'cerebellum neurons' are in `subcluster_cell_type-Di- and mesencephalon excitatory neurons.ipynb` and `subcluster_cell_type-cerebellumneuron.ipynb`. Change main cell type input for susbclustering of other main cell types.


### Plot UMAP and cell type map
- [UMAP](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.umap.html#scanpy-tl-umap) is plotted colored by subclusters in `plot_umap_cellmap.ipynb`.
- Spatial cell-type maps are plotted colored by subclusters of 20 mouse CNS slices. Each dot represents one cell. 
