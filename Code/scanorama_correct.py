'''
This script inputs a directory and a list of single-cell count matrices in the directory and integrates and batch corrects the datasets using Scanorama.  

*All datasets will be filtered for a minimum of 200 genes per cell and a minimum of 30 cells per gene.  Certain other quality control steps are skipped because the datasets we are using were preprocessed.  All except the Lawlor et al. data set had mitochondrial reads and cells with too few reads or high % of mitochondrial reads removed.  The Lawlor et al data was scrubbed for low read counts but not explicitly for mitochondrial read percentages.  The authors commented that the cells with low read counts also had some of the highest mitochondrial percentages. 

**Some datasets had been previously normalized or log transformed which will impact the efficacy of scVI as it requires raw counts.  For more details on each data set please refer to Baron et al. (2016), Muraro et al. (2016), Grun et al. (2016), Lawlor er al. (2017), and Segerstolpe et al. (2016).  
'''

#######################################################High-Level Variables##################################################################

#Set directories
raw_dir = '/home/conerym/BSTA_787/raw_pancreas_data/'
process_dir = '/home/conerym/BSTA_787/processed_pancreas_data/scanorama/'

#Set list of data sets
data_names = [
    raw_dir + 'pancreas_inDrop.txt',
    raw_dir + 'pancreas_multi_celseq2_expression_matrix.txt',
    raw_dir + 'pancreas_multi_celseq_expression_matrix.txt',
    raw_dir + 'pancreas_multi_fluidigmc1_expression_matrix.txt',
    raw_dir + 'pancreas_multi_smartseq2_expression_matrix.txt',
]

#Set desired number of clusters
cluster_count = 13 #Determined number of desired clusters for pancreas
#alpha, beta, gamma, delta, epsilon, acinar, activated_stellate, ductal, endothelial, macrophage, mast, quiescent_stellate, and schwann

#Read in cell type annotations
TYPE_annotations = '/home/conerym/BSTA_787/raw_pancreas_data/cell_type_annotations.tsv'
#Set marker genes
marker_genes_dict={'Alpha': ['GCG'],'Beta': ['INS'],'Gamma': ['PPY'], 'Delta': ['SST'], 'Ductal': ['KRT19'], 'Acinar': ['PRSS1'], 'Stellate': ['COL1A1']}

#Make results files
scanorama_file = process_dir + "scanorama_corrected.h5ad"

########################################################Import Modules######################################################################################


#import modules that will be needed
import sys
from PIL import Image
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import scanpy.external as sce
import matplotlib.pyplot as plt
from copy import deepcopy
from shutil import move
import warnings

'''Machine learning and single-cell packages'''
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
import scanpy as sc
from anndata import AnnData

'''Merging pacakges'''
import scanorama

############################################################Define Useful Functions#####################################################################

def find_resolution(adata_, n_clusters, random = 0): 
    adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        sc.tl.louvain(adata, resolution = current_res, random_state = random)
        labels = adata.obs['louvain']
        obtained_clusters = len(np.unique(labels))
        
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        
        iteration = iteration + 1
        
    return current_res

############################################################### Scanorama ##################################################################

'''We follow the instruction used in the original Hie et al. (2019) paper for l2-normalization'''

#Read in data files and store as a list of Anndata objects.  Execute needed qc on the objects.  
#Also make a list of shared genes across all cells
datasets = []
for i in range(len(data_names)):
    #Read in datasets and remove rarely expressed genes and cells with too few genes
    datasets.append(sc.read(data_names[i], delimiter="\t", first_column_names=True))
    datasets[i]=datasets[i].transpose()
    sc.pp.filter_cells(datasets[i], min_genes=200)
    sc.pp.filter_genes(datasets[i], min_cells=30)
    
    #l2 Normalize data
    datasets[i].X = sp.sparse.csr_matrix(datasets[i].X)
    sce.tl.sam(datasets[i], k=20, distance='euclidean', projection='None', standardization='Normalizer', inplace = True)
    datasets[i].X = sp.sparse.csr_matrix.toarray(datasets[i].X)
    #sc.pp.normalize_total(datasets[i], target_sum=1e4)
        
    #Make gene lists
    if i != 0:
        genes = [value for value in genes if value in list(datasets[i].obs_names)]
    else:
        genes = list(datasets[i].obs_names)

#Scanorama data integration and batch correction
scanorama_corrected = scanorama.correct_scanpy(datasets, return_dimred=True)

#Create adata by concatenation (This step is not readily expandable to a greater/smaller number of datasets
adata = scanorama_corrected[0].concatenate(scanorama_corrected[1], scanorama_corrected[2], scanorama_corrected[3], scanorama_corrected[4], batch_categories=["inDrop", "celseq2", "celseq", "fluidigmc1", "smartseq2"])
print(adata)
type(adata)

#Read in annotations
type_annotations = {}
a_file = open(TYPE_annotations, encoding='utf-8')
for line in a_file:
    key, value = line.split()
    type_annotations[key] = value
del line
del key
del value
#Add a new T2D status column using our t2d annotations file
adata.obs['cell type'] = adata.obs_names.map(type_annotations).astype('category')

#Execute PCA on batch corrected data
sc.tl.pca(adata, svd_solver='arpack')
#Make PCA plot and elbow plot
sc.pl.pca(adata, color='batch', save='_scanorama.png')
sc.pl.pca_variance_ratio(adata, save='_scanorama.png')

#Calculate neighbors 
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
#Cluster using Louvain Algorithm
res_scanorama = find_resolution(adata, cluster_count)
sc.tl.louvain(adata, resolution=res_scanorama)
adata.obs['cluster assignment'] = adata.obs['louvain']
print(adata)
#Make UMAP plot
sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch', 'cell type', 'cluster assignment'], show=False, ncols=1, save='_scanorama.png') #Colored by batch and cluster

#Make dotplot
sc.pl.dotplot(adata, marker_genes_dict, 'cluster assignment', dendrogram=True, save= 'scanorama.png')

#save to file
adata.write_csvs(process_dir)
adata.write_h5ad(scanorama_file)




