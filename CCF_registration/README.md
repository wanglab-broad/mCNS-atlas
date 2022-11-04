### CCF registration
- We performed registration of each STARmap PLUS slice with [Allen CCF v3](https://pubmed.ncbi.nlm.nih.gov/32386544/)
- To match each STARmap PLUS slice to their corresponding CCF slices, we employed the [AP_histology package](https://github.com/petersaj/AP_histology) and made some adjustments.
- The analysis included the following steps:
  -   [Annotate](AP_histology-adapted-registration/read_allin_ccf_mapping_pipeline.m)
  -   [Display](AP_histology-adapted-registration/retrieve_previous_mapping_multiple.m)
  -   [Get cell region ID labels](AP_histology-adapted-registration/get_histology_annotation_on_raw.m)
  -   [Get cell region parent ID labels](AP_histology-adapted-registration/readin_ccf_annotation.mlx)
