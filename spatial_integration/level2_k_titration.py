import scanpy as sc
import numpy as np
import spin_cat
import argparse


def main(num_nbrs):

    print('Loading data', flush=True)
    adata = sc.read_h5ad('/stanley/WangLab/kamal/data/integrated/spatial/atlas/atlas_spin_concat_nocombat_res015.h5ad')
    samples = adata.obs['sample'].unique()
    xkey = 'row'
    ykey = 'col'

    ress = [
        1.7, # 0 Ctx 1 (separate L4 latent branches, which might be middle vs ends)
        1.7, # 1 Midbrain 1 (separate OB-adjacent and SCs vs SCm)
        0.4, # 2 Cb 1
        0.4, # 3 Str (3 levels, D to V, and ideally fiber streaks)
        0.1, # 4 Ob 1
        0.2, # 5 fiber tracts (dense region of neurogenesis in wm by thal in sagittal3)
        0.2, # 6 meninges (separate out HPF)
        0.4, # 7 Ctx 2 (Amygldala, Piriform)
        0.1, # 8 Ob 2
        0.5, # 9 Thal
        0.3, # 10 Cb 2
        0.6, # 11 choroid (separate out MHb)
        0.3, # 12 DG (two very separate clusters, so will try and segment it)
        0.1, # 13 CA
        0.05, # 14 Midbrain 2 (just check marker genes to justify the resulting 3 distinct clusters)
        0.0, # 15 spinal cord (Dorsal)
        0.3, # 16 Ctx 3 (should have layers)
    ]

    topregion_key = 'region'
    adata.obs['topregion'] = adata.obs[topregion_key]
    subregion_key = 'subregion'
    subpca_key = 'X_pca_spin_sub'
    subumap_key = 'X_umap_spin_sub'
    subdatas = []
    final_region_key = 'region'
    count_thresh = 100

    # Filter small regions (selected by hand)
    print('Filtering out small regions', flush=True)
    good_regions = [str(region) for region in range(17)]
    adata = adata[np.isin(adata.obs[topregion_key], good_regions)]
    print(adata)
    print(adata.obs[topregion_key])
    regions = adata.obs[topregion_key].unique()
    print('num_nbrs = ', num_nbrs, flush=True)

    print('Subclustering', flush=True)
    for i in range(len(regions)):
        
        # Select region
        region = regions[i]
        print(region, flush=True)
        subdata = adata[adata.obs[topregion_key]==region]
        
        # Split into separate samples for spatial integration and filter
        subdata_list = [subdata[subdata.obs['sample']==sample] for sample in samples]
        
        # Filter subdatas
        subdata_list = [subdata_list[i] for i in range(len(subdata_list)) if len(subdata_list[i])>=count_thresh] # filter based on counts for each region (30 minimum for SPIN)
        sample_list = [samples[i] for i in range(len(subdata_list)) if len(subdata_list[i])>=count_thresh]
        if len(sample_list) == 0: # skip region if no sample has >=count_thresh of that region in it
            continue

        # Cluster the given region across samples
        pooldata, _ = spin_cat.integrate(subdata_list, sample_list, pca_key=subpca_key, xkey=xkey, ykey=ykey, svd_solver='randomized',
                                        num_nbrs=num_nbrs)
        spin_cat.cluster(pooldata, pca_key=subpca_key, region_key=subregion_key, umap_key=subumap_key, res=ress[i])

        # Accumulate results
        subdatas.append(pooldata)

        # Write results to recluster later (integration only needed once)
        print('Writing to disk', flush=True)
        pooldata.write(f'/stanley/WangLab/kamal/data/integrated/spatial/atlas/nbr_titration/atlas_spin_concat_nocombat_subclustered_topregion{region}_{num_nbrs}nbrs.h5ad')

    # Combine subclustered data
    adata_sub = sc.concat(subdatas)

    # Create final region labels
    adata_sub.obs[final_region_key] = [adata_sub.obs[topregion_key][i]+'_'+adata_sub.obs[subregion_key][i] for i in range(len(adata_sub))]
    adata_sub.obs[final_region_key] = adata_sub.obs[final_region_key].astype('category')

    # Write data to disk
    print(f'Writing data to /stanley/WangLab/kamal/data/integrated/spatial/atlas/nbr_titration/atlas_spin_concat_nocombat_subclustered_{num_nbrs}nbrs.h5ad', flush=True)
    adata_sub.write(f'/stanley/WangLab/kamal/data/integrated/spatial/atlas/nbr_titration/atlas_spin_concat_nocombat_subclustered_{num_nbrs}nbrs.h5ad')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--num_nbrs', type=int)
    args = parser.parse_args()

    main(args.num_nbrs)