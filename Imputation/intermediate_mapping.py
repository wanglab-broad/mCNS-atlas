#!/usr/bin/env python
# coding: utf-8

import sys, os
tile = int(sys.argv[1])
print(f'tile:{tile}')
import scanpy.external as sce
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
from scipy.stats import spearmanr,pearsonr
from sklearn.preprocessing import MinMaxScaler,Normalizer,StandardScaler,quantile_transform
import os
print(f'tile:{tile}')
### choose imputated gene
genelist=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/imputation/input_data/genelist.csv',index_col=None,header=None)
predict_genelist=genelist.iloc[tile-1][0]
print(f'tile:{tile}')
###read STARmap data
adata=sc.read_h5ad('/n/holystore01/LABS/jialiu_lab/Lab/imputation/input_data/adata_afterbatch.h5ad')
print(f'tile:{tile}')
adata_predict=adata[:,[x not in predict_genelist for x  in adata.var.index]]
adata_org=adata.copy()
adata=adata_predict
###make gene name CAPTION
adata.var.index=adata.var.index.str.upper()


### read single cell data
adata_sc=sc.read_h5ad('/n/holystore01/LABS/jialiu_lab/Lab/imputation/input_data/scrna.h5ad')
adata_sc.var.index=adata_sc.var.index.str.upper()
adata_sc.var_names_make_unique()

list_of_variable_names=adata.var.index.intersection(adata_sc.var.index)
adata_subset = adata[:, list_of_variable_names]

list_of_variable_names=adata_sc.var.index.intersection(adata.var.index)
adata_sc_subset = adata_sc[:, list_of_variable_names]



################ harmony
# Normalization scaling sc
sc.pp.normalize_total(adata_sc_subset)
sc.pp.log1p(adata_sc_subset)

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata_sc_subset, percent_top=None, inplace=True)

# Scale data to unit variance and zero mean
sc.pp.regress_out(adata_sc_subset, ['total_counts'])
sc.pp.scale(adata_sc_subset)

subset=adata_sc_subset.obs.drop(adata_sc_subset[adata_sc_subset.obs['TaxonomyRank2']=='PNS neurons'].obs.index)
adata_sc_subset=adata_sc_subset[subset.index,:]

subset=adata_sc_subset.obs.drop(adata_sc_subset[adata_sc_subset.obs['Description']=='Schwann cells'].obs.index)
adata_sc_subset=adata_sc_subset[subset.index,:]

subset=adata_sc_subset.obs.drop(adata_sc_subset[adata_sc_subset.obs['Description']=='Satellite glia'].obs.index)
adata_sc_subset=adata_sc_subset[subset.index,:]

subset=adata_sc_subset.obs.drop(adata_sc_subset[adata_sc_subset.obs['Description']=='Enteric glia'].obs.index)
adata_sc_subset=adata_sc_subset[subset.index,:]

combine_adata=adata_subset.concatenate(adata_sc_subset,batch_key='dataset', batch_categories=['st','scrna'])

sc.tl.pca(combine_adata,n_comps=100)
sce.pp.harmony_integrate(combine_adata, 'dataset')
print('harmony finished')
sc.pp.neighbors(combine_adata,n_neighbors=50,use_rep='X_pca_harmony')

sc.tl.umap(combine_adata)



### save umap plot
cluster_pl = sns.color_palette("tab20",4)
cluster_pl[0]=np.array([252,40,3])/255
cluster_pl[1]=np.array([3, 169, 252])/255
cluster_pl[2]=np.array([4,217,61])/255
sns.palplot(cluster_pl)
fig,ax = plt.subplots(figsize=(8,8))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(7)
ax=sc.pl.umap(combine_adata, color='dataset',legend_fontsize=13,ax=ax,show=False)
ax.set_title('Harmony')
ax.title.set_fontsize(20)
fig.savefig('/n/holystore01/LABS/jialiu_lab/Lab/imputation/output/figure/umap_'+predict_genelist+'.pdf')


### save combined adata
# combine_adata.write_h5ad('/n/holystore01/LABS/jialiu_lab/Lab/imputation/output/data/adata_harmony_integrate_'+predict_genelist+'.h5ad')


########prepare for imputation
combine_adata_st=combine_adata[combine_adata.obs['dataset']=='st',:].copy()
combine_adata_sc=combine_adata[combine_adata.obs['dataset']=='scrna',:].copy()

adata_sc_imputation=adata_sc.copy()

subset=adata_sc_imputation.obs.drop(adata_sc_imputation[adata_sc_imputation.obs['TaxonomyRank2']=='PNS neurons'].obs.index)
adata_sc_imputation=adata_sc_imputation[subset.index,:]
subset=adata_sc_imputation.obs.drop(adata_sc_imputation[adata_sc_imputation.obs['Description']=='Schwann cells'].obs.index)
adata_sc_imputation=adata_sc_imputation[subset.index,:]
subset=adata_sc_imputation.obs.drop(adata_sc_imputation[adata_sc_imputation.obs['Description']=='Satellite glia'].obs.index)
adata_sc_imputation=adata_sc_imputation[subset.index,:]
subset=adata_sc_imputation.obs.drop(adata_sc_imputation[adata_sc_imputation.obs['Description']=='Enteric glia'].obs.index)
adata_sc_imputation=adata_sc_imputation[subset.index,:]

adata_sc_imputation_subset=adata_sc_imputation[:,[list(adata_sc_imputation.var.index).index(predict_genelist)]]
singlecell_rawexpr=adata_sc_imputation_subset.X
shape1=singlecell_rawexpr.shape[1]

adata=sc.read_h5ad('/n/holystore01/LABS/jialiu_lab/Lab/imputation/input_data/adata_afterbatch.h5ad')
adata_real=adata[:,[list(adata.var.index).index(predict_genelist)]]
adata_realexpr=adata_real.X


pearsonrlist={'gene':[],'value':[],'pvalue':[],'n_neighbors':[]}

for n_neighbors in tqdm([5,10,20,40,50,100,125,200,250,500,1000,2000]):
    neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(combine_adata_sc.obsm['X_umap'])
    distances, indices = neigh.kneighbors(combine_adata_st.obsm['X_umap'])

    indices_all=np.concatenate(indices)

    adata_st_imputation=[]
    del adata_st_imputation

    shape1=singlecell_rawexpr.shape[1]
    for i in range(0,indices_all.shape[0],100000):

        test1=singlecell_rawexpr[indices_all[i:i+100000]].toarray()
        test2=test1.reshape(int(test1.shape[0]/n_neighbors),n_neighbors,shape1)
        test3=np.mean(test2,axis=1)
        if i==0:
            adata_st_imputation=test3
        else:
            adata_st_imputation=np.concatenate((adata_st_imputation,test3),axis=0)

    adata_imput_expr=adata_st_imputation



    focus_gene=0
    gene_name=predict_genelist
    dis1=adata_imput_expr
    dis2=np.array(adata_realexpr[range(adata_subset.shape[0]),:]).copy()
    pearsonrlist['gene'].append(gene_name)
    pearsonrlist['value'].append(pearsonr(dis1.flatten(),dis2.flatten())[0])
    pearsonrlist['pvalue'].append(pearsonr(dis1.flatten(),dis2.flatten())[1])
    pearsonrlist['n_neighbors'].append(n_neighbors)
pd.DataFrame(pearsonrlist).to_csv('/n/holystore01/LABS/jialiu_lab/Lab/imputation/output/data/pearsonrlist'+predict_genelist+'.csv')