import numpy as np 
import pandas as pd 
import scanpy as sc
import sys 
from statsmodels import robust 
import sys 
from scipy import sparse, io
import matplotlib.pyplot as plt 
import os.path
import scipy as sci
import anndata
import bbknn
import anndata

os.chdir('/home/ngr18/covid/sanger_data/sanger_data/built_data')
sanger = sc.read_h5ad('sanger.h5ad')

os.chdir('/home/ngr18/covid/sanger_data/cambridge_data/built_data')
cam = sc.read_h5ad('cam.h5ad')

os.chdir('/home/ngr18/covid/sanger_data/ncl_data')
ncl = sc.read_h5ad('ncl.h5ad')

total = anndata.concat([cam, sanger, ncl], index_unique = None)

total.var['feature_types'] = np.where(total.var.index.str.contains('AB_'), 'Antibody Capture',
                                      'Gene Expression')

total.obs.rename(columns = {"patient_id": "sample_id"}, inplace = True)

os.chdir('/home/ngr18/covid/sanger_data/')
total.write('combined_dec.h5ad')

rna = total[:, total.var["feature_types"] == "Gene Expression"].copy()

sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)

rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

rna = rna[rna.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(rna, target_sum=1e4)

sc.pp.log1p(rna)

sc.pp.highly_variable_genes(rna, n_top_genes = 3000, flavor = 'seurat')

rna.raw = rna

rna = rna[:, rna.var.highly_variable]

sc.pp.scale(rna, max_value=10)

sc.tl.pca(rna, svd_solver='arpack')

sc.external.pp.harmony_integrate(rna, 'sample_id', basis='X_pca', adjusted_basis='X_pca_harmony')

sc.pp.neighbors(rna, n_neighbors=10, n_pcs=40, use_rep = 'X_pca_harmony')

sc.tl.umap(rna)

sc.tl.leiden(rna, resolution = 3)

os.chdir('/home/ngr18/covid/sanger_data/')
rna.write('combined_harmonyv3.h5ad')


