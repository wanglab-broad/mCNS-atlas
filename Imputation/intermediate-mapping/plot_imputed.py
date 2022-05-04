#!/usr/bin/env python
# coding: utf-8

import sys, os
tile = int(sys.argv[1])

print(f'tile:{tile}')

import os
import warnings 
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from anndata import AnnData
from pathlib import Path

from pathlib import Path


combine_adata_st_obs=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/combine_adata_st_obs.csv',index_col=0)

celltypeannotaion='/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data_save/'
for wellname in combine_adata_st_obs['sample'].unique():
    Path(celltypeannotaion+wellname).mkdir(parents=True, exist_ok=True)


good_gene=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/good_gene.csv')


well_name=combine_adata_st_obs['sample'].unique()[tile]
print(well_name)
loaded = np.load('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data_save/'+well_name+'_imputed_score.npz')
st_imputation_well=loaded['a']

combine_adata_st_sub_obs=combine_adata_st_obs.loc[combine_adata_st_obs['sample']==well_name,:].copy()
adata_imputed=AnnData(st_imputation_well,var=pd.DataFrame(good_gene['Gene']),obs=combine_adata_st_sub_obs)

# Normalization scaling
print('normalize'+well_name)
sc.pp.normalize_total(adata_imputed)
sc.pp.log1p(adata_imputed)

# Scale data to unit variance and zero mean
sc.pp.regress_out(adata_imputed, ['total_counts'])
sc.pp.scale(adata_imputed)


for idx in range(11844):
    plt.figure(figsize=(max(adata_imputed.obs['col'])/1000,max(adata_imputed.obs['row'])/1000))
    plt.scatter(adata_imputed.obs['col'],adata_imputed.obs['row'],
               s=6,
                vmin=-5,vmax=5,
               c=adata_imputed.X[:,idx],
               cmap='coolwarm')
    plt.savefig(celltypeannotaion+well_name+'/'+adata_imputed.var['Gene'][idx]+'.png',
               dpi=300)



