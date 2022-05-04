#!/usr/bin/env python
# coding: utf-8

import sys, os
tile = int(sys.argv[1])

print(f'tile:{tile}')

import os
import warnings 
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd


distances = pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/distances.csv',index_col=0)
indices = pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/indices.csv',index_col=0)
adata_sc_imputation_X = pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/adata_sc_imputation_X.csv',index_col=0)
combine_adata_st_obs = pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/combine_adata_st_obs.csv')
distances=np.array(distances)
indices=np.array(indices)
adata_sc_imputation_X=np.array(adata_sc_imputation_X)

# distances=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/distances.csv')
# indices=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/indices.csv')
# adata_sc_imputation_X=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/adata_sc_imputation_X.csv')
# combine_adata_st_obs=pd.read_csv('/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data/combine_adata_st_obs.csv.csv')
specific_batch=np.where(combine_adata_st_obs['sample']==combine_adata_st_obs['sample'].unique()[tile])


st_imputation_well=np.zeros((specific_batch[0].shape[0],11844))
for indi,i in enumerate(specific_batch[0]):
    test=1/distances[i]#
    test=test/test.sum()
    test1=np.dot(test, adata_sc_imputation_X[indices[i],:])
    st_imputation_well[indi,:]=np.array(test1)

well_name=combine_adata_st_obs['sample'].unique()[tile]
np.savez_compressed(f'/n/holystore01/LABS/jialiu_lab/Lab/final_impu/data_save/{well_name}_imputed_score.npz',
                    a=st_imputation_well)

