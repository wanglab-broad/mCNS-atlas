Scripts
spin_cat.py: core spatial integration script used to spatially cluster across multiple slices
level2_k_titration.py: spatial subclustering script
level2_k_titration.sh: shell script for submitting parallel runs of level2_k_titration.py

Notebooks
level2_param_selection.ipynb: for manual selection of subclustering parameters, such as neighborhood sizes k and Leiden resolutions
region1_level3_clustering.ipynb: for manual selection (and execution) of subsubclustering parameters, such as neighborhood sizes k and Leiden resolutions
region5_level3_clustering.ipynb: for manual selection (and execution) of subsubclustering parameters, such as neighborhood sizes k and Leiden resolutions
region11_level3_clustering.ipynb: for manual selection (and execution) of subsubclustering parameters, such as neighborhood sizes k and Leiden resolutions
combine_subclustered_adatas.ipynb: for combining clustering results across levels 1-3 into a single .h5ad file
