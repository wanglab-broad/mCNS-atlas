#TODO: try finding hvgs AFTER concatting

import scanpy as sc
import numpy as np
import anndata as ad
import squidpy as sq
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import argparse


def spatially_variable_genes(adata, method='moran', num_genes=1000):

    if method=='moran':
        sq.gr.spatial_neighbors(adata, key_added='spatial')
        sq.gr.spatial_autocorr(adata, mode='moran')
        adata.var = pd.concat([adata.var, adata.uns['moranI']], axis=1)
        thresh = adata.var['I'].sort_values(ascending=False)[num_genes]
        adata.var['spatially_variable'] = adata.var['I'] > thresh # possible to get more than num_genes


def smooth(adata, num_nbrs=30, num_samples=10, xkey='x', ykey='y', num_drop=None,
           subsample=False, use_highly_variable=False, use_spatially_variable=False):

    # Find nbrhds
    coordinates = adata.obs[[xkey, ykey]]
    nbrs = NearestNeighbors(n_neighbors=num_nbrs)
    nbrs.fit(coordinates)
    _, nbr_idxs = nbrs.kneighbors(coordinates)

    # Subsample nbrhds
    if num_drop:
        nbr_idxs_sampled = nbr_idxs[:,num_drop:]
    elif subsample:
        nbr_idxs_sampled = [np.random.choice(idxs, size=num_samples) for idxs in nbr_idxs]
    else:
        nbr_idxs_sampled = nbr_idxs

    # Isolate highly/spatially variable genes
    if use_highly_variable:
        X = adata.X[:, adata.var['highly_variable']]
    elif use_spatially_variable:
        X = adata.X[:, adata.var['spatially_variable']]
    else:
        X = adata.X

    # Concatenate nbr vectors
    X_smooth = []
    for nth_nbrs in np.array(nbr_idxs_sampled).T:
        X_smooth.append(X[nth_nbrs])

    return np.hstack(X_smooth)


def integrate(adatas, batch_labels, num_nbrs=30, num_samples=10, xkey='x', ykey='y',
              num_pcs=50, svd_solver='randomized', batch_key='batch', num_drop=None,
              pca_key='X_pca_spin', use_highly_variable=False, num_hvgs=None,
              subsample=False, use_spatially_variable=False, svg_method='moran',
              num_genes=1000):

    # Keep shared genes
    gene_sets = [set(adata.var_names) for adata in adatas]
    shared_genes = list(set.intersection(*gene_sets))
    adatas = [adata[:,shared_genes] for adata in adatas]

    # Concatenate datasets
    pooldata = ad.concat(adatas, keys=batch_labels, label=batch_key, join='outer')

    # # Batch correct
    # if len(adatas) > 1:
    #     print('Batch correcting', flush=True)
    #     sc.pp.combat(pooldata, key=batch_key)

    # Find variable genes
    if use_highly_variable:
        print('Calculating HVGs', flush=True)
        print(f'num_hvgs: {num_hvgs}', flush=True)
        sc.pp.highly_variable_genes(pooldata, batch_key=batch_key, n_top_genes=num_hvgs)
    elif use_spatially_variable:
        print('Calculating SVGs', flush=True)
        spatially_variable_genes(pooldata, method=svg_method, num_genes=num_genes)

    # Spatially smooth
    X_smooth = []
    for i in range(len(adatas)):
        subdata = pooldata[pooldata.obs[batch_key]==batch_labels[i]]
        X_smooth.append(
            smooth(
                subdata,
                num_nbrs=num_nbrs,
                num_samples=num_samples,
                xkey=xkey,
                ykey=ykey,
                use_highly_variable=use_highly_variable,
                use_spatially_variable=use_spatially_variable,
                num_drop=num_drop,
                subsample=subsample
            )
        )
    X_smooth = np.vstack(X_smooth)

    # Run PCA
    pca = PCA(n_components=num_pcs, svd_solver=svd_solver)
    pooldata.obsm[pca_key] = pca.fit_transform(X_smooth)

    # Integrate PCs
    if len(adatas) > 1:
        sc.external.pp.harmony_integrate(
            pooldata,
            batch_key,
            basis=pca_key,
            adjusted_basis=pca_key
        )

    return pooldata, pca


def cluster(pooldata, pca_key='X_pca_spin', region_key='region', umap_key='X_umap_spin',
            res=0.5):

    sc.pp.neighbors(pooldata, use_rep=pca_key, key_added=region_key)
    umap = sc.tl.umap(pooldata, neighbors_key=region_key, copy=True).obsm['X_umap']
    pooldata.obsm[umap_key] = umap
    sc.tl.leiden(pooldata, resolution=res, key_added=region_key, neighbors_key=region_key)

    return pooldata


def main(adata_paths, batch_keys, num_nbrs, num_samples, xkey, ykey, res, write_path,
         num_pcs, svd_solver, batch_key, region_key, pca_key, umap_key, num_drop,
         use_hvgs, num_hvgs):

    # Import data
    adatas = [sc.read_h5ad(path) for path in adata_paths]

    # Integrate spatial features across samples
    print('Integrating spatial features\n', flush=True)
    pooldata, _ = integrate(
        adatas,
        batch_keys,
        num_nbrs=num_nbrs,
        num_samples=num_samples,
        xkey=xkey,
        ykey=ykey,
        num_pcs=num_pcs,
        svd_solver=svd_solver,
        batch_key=batch_key,
        pca_key=pca_key,
        num_drop=num_drop,
        use_highly_variable=use_hvgs,
        num_hvgs=num_hvgs
    )

    # Cluster integrated samples to find regions
    print('\nClustering on spatial features\n', flush=True)
    pooldata = cluster(pooldata, pca_key=pca_key, region_key=region_key,
                       umap_key=umap_key, res=res)

    # Write output
    pooldata.write(write_path)
    print(f'\tOutput written to {write_path}', flush=True)


if __name__=='__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='SPatially INtegrate and cluster one or more spatial omics samples'
    )
    parser.add_argument('--adata_paths', '-p', nargs='+', help='Path to each dataset')
    parser.add_argument('--batch_keys', '-b', nargs='+', help='List of dataset labels')
    parser.add_argument('--num_nbrs', '-k', type=int, default=30)
    parser.add_argument('--num_samples', '-s', type=int, default=10)
    parser.add_argument('--xkey', '-x', type=str, default='x')
    parser.add_argument('--ykey', '-y', type=str, default='y')
    parser.add_argument('--res', type=float, default=0.5)
    parser.add_argument('--write_path', '-w', help='Path to output file')
    parser.add_argument('--num_pcs', type=int, default=50)
    parser.add_argument('--svd_solver', default='randomized')
    parser.add_argument('--batch_key', default='batch_key')
    parser.add_argument('--region_key', default='region')
    parser.add_argument('--pca_key', default='X_pca_spin')
    parser.add_argument('--umap_key', default='X_umap_spin')
    parser.add_argument('--num_drop', default=None)
    parser.add_argument('--use_hvgs', action='store_true')
    parser.add_argument('--num_hvgs', type=int, default=None)
    args = parser.parse_args()

    # Run integration with parsed arguments
    main(args.adata_paths, args.batch_keys, args.num_nbrs, args.num_samples, args.xkey,
         args.ykey, args.res, args.write_path, args.num_pcs, args.svd_solver,
         args.batch_key, args.region_key, args.pca_key, args.umap_key, args.num_drop,
         args.use_hvgs, args.num_hvgs)
