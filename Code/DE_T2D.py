'''
This script tests for DE genes within each cluster between known T2D cases and controls.  
'''

#######################################################Testing Variables##################################################################
'''

#Set input files and output_dir
T2D_annotations = '/home/conerym/BSTA_787/raw_pancreas_data/diabetes_annotations.tsv'
output_dir= '/home/conerym/BSTA_787/processed_pancreas_data/scanorama/'
input_file = output_dir + "scanorama_corrected.h5ad"

'''
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
import random
import csv
import statsmodels.stats.multitest as multi

'''Machine learning and single-cell packages'''
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
import scanpy as sc
from anndata import AnnData

#########################################################Main Thread####################################

def main_thread(output_dir, input_file, T2D_annotations):
    
    #Read in data file
    adata = sc.read_h5ad(input_file)
    
    #Read in annotations
    t2d_annotations = {}
    a_file = open(T2D_annotations, encoding='utf-8')
    for line in a_file:
        key, value = line.split()
        t2d_annotations[key] = value
    del line
    del key
    del value
    #Add a new T2D status column using our t2d annotations file
    adata.obs['T2D Status'] = adata.obs_names.map(t2d_annotations).astype('category')
    
    ####################### Cell Types #################################
    
    #Make output file
    if output_dir[-1] == "/":
        output_loc = output_dir + "T2D_DE_cell_type_results.tsv"
    else:
        output_loc = output_dir + "/" + "T2D_DE_cell_type_results.tsv" 
    #Make list to hold output for export
    genes = list(set(adata.var_names.tolist()))
    genes.insert(0, 'gene')
    #convert every entry to a list to facilitate appending
    for i in range(len(genes)):
        genes[i] = [genes[i]]

    #Cycle through all clusters and conduct the DE analyses between cases/controls
    for cell_type in list(set(adata.obs["cell type"].tolist())):  #This avoids the nan type
        
        #Create two arrays by subsetting for T2D status and cluster
        t2d_adata = adata[adata.obs['T2D Status'] == '1']
        t2d_adata = t2d_adata[t2d_adata.obs['cell type'] == cell_type]
        t2d_adata.obs["value"] = 0
        t2d_array = sp.sparse.csr_matrix.toarray(t2d_adata.X)
        nd_adata = adata[adata.obs['T2D Status'] == '0']
        nd_adata = nd_adata[nd_adata.obs['cell type'] == cell_type]
        nd_adata.obs["value"] = 0
        nd_array = sp.sparse.csr_matrix.toarray(nd_adata.X)
        
        #Calculate p-values for each gene 
        DE_results = execute_wilcox(nd_array, t2d_array)
        #Insert headers onto DE_results
        DE_results.insert(0, [str(cell_type) + '_stat', str(cell_type) + '_pval'])
        #Cycle through DE results and extend entries to genes list
        for i in range(len(DE_results)):
            genes[i].extend([DE_results[i][0], DE_results[i][1]])

    #write file
    with open(output_loc, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in genes:
            writer.writerow(i)
        writeFile.close()
    
    ####################### Clusters #################################
    
    #Make output file
    if output_dir[-1] == "/":
        output_loc = output_dir + "T2D_DE_cluster_results.tsv"
    else:
        output_loc = output_dir + "/" + "T2D_DE_cluster_results.tsv" 
    #Make list to hold output for export
    genes = list(set(adata.var_names.tolist()))
    genes.insert(0, 'gene')
    #convert every entry to a list to facilitate appending
    for i in range(len(genes)):
        genes[i] = [genes[i]]

    #Cycle through all clusters and conduct the DE analyses between cases/controls
    for cluster in list(set(adata.obs["cluster assignment"].tolist())):
        
        #Create two arrays by subsetting for T2D status and cluster
        t2d_adata = adata[adata.obs['T2D Status'] == '1']
        t2d_adata = t2d_adata[t2d_adata.obs['cluster assignment'] == cluster]
        t2d_adata.obs["value"] = 0
        t2d_array = sp.sparse.csr_matrix.toarray(t2d_adata.X)
        nd_adata = adata[adata.obs['T2D Status'] == '0']
        nd_adata = nd_adata[nd_adata.obs['cluster assignment'] == cluster]
        nd_adata.obs["value"] = 0
        nd_array = sp.sparse.csr_matrix.toarray(nd_adata.X)
        
        #Calculate p-values for each gene 
        DE_results = execute_wilcox(nd_array, t2d_array)
        #Insert headers onto DE_results
        DE_results.insert(0, [cluster + '_stat', cluster + '_pval'])
        #Cycle through DE results and extend entries to genes list
        for i in range(len(DE_results)):
            genes[i].extend([DE_results[i][0], DE_results[i][1]])

    #write file
    with open(output_loc, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in genes:
            writer.writerow(i)
        writeFile.close()

#####################################################Statistics Functions############################################

'''
This function executes the Wilcoxon tests it takes in two arrays that are assumed to have genes along the columns and cells down the rows.  It returns a list of p-values
'''

def execute_wilcox(array_1, array_2):
    
    #Make a blank list to store the p-values
    p_values = []
    
    #Iterate over the columns (genes) in the arrays (it should be the same in both!)
    for i in range(array_1.shape[1]):
        #Execute Wilcoxon test and add p_value to p_values
        p_values.append(sp.stats.ranksums(array_1.T[i], array_2.T[i]))
    #Return p_values
    return p_values

#########################################################Read in Arguments#############################################

if len(sys.argv) == 4:         
    #Call main thread
    main_thread(output_dir=sys.argv[1], input_file=sys.argv[2], T2D_annotations=sys.argv[3])
else:
    print("Incorrect number of parameters given.  Please try again.")
    print("Paramters should be the following: output directory, input file(full path), T2D annotations (fullpath)")

