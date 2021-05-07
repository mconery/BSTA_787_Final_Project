'''
This script inputs a h5ad file with 'cluster assignment' and 'batch' annotations, subsets the data into two groups, and then uses a Mann-Whitney-Wilcoxon Rank-Sum test to test for the number of DE genes between the two groups.  It executes this process twice, once randomly and once in an ordered manner, for an input number of iterations.  Both the random and ordered methods for subsetting the data ensure that an approximately even number of cells from each cluster will be allocated to each subset.  The random allocation within each cluster is truly random, so that each subset will wind up with an even 50-50 distribution of each batch.  The ordered allocation, which assumes the input AnnData object is presorted in batch order, picks a random starting point within each clustered set of cells and consecutively allocates half to one subset and the other to the second subset.  This method ensures crude stratification of each subset by batch.  

Paramters should be the following: naming convention, output directory, input file(full path), sample size (optional/default=100), alpha (optional/default=0.05).  If either optional parameter is given both must be given.  
'''

#######################################################Testing Variables##################################################################
'''

#Set directories
process_dir = '/home/conerym/BSTA_787/processed_pancreas_data/scanorama/'
#Set naming convention and alpha
naming_conv = 'scanorama'
alpha=0.05
sample_size=100
#Set input file and output_dir
output_dir=process_dir
input_file = process_dir + "scanorama_corrected.h5ad"

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

def main_thread(naming_conv, output_dir, input_file, sample_size, alpha):
    
    #Read in data file
    adata = sc.read_h5ad(input_file)
    
    #Make output files
    if output_dir[-1] == "/":
        ordered_file = output_dir + naming_conv + "_ordered_pvals.tsv"
        random_file =  output_dir + naming_conv + "_random_pvals.tsv"
    else:
        ordered_file = output_dir + "/" + naming_conv + "_ordered_pvals.tsv"
        random_file =  output_dir + "/" + naming_conv + "_random_pvals.tsv"    
    
    #Create blank lists to hold p_value lists
    ordered_p_vals=[]
    random_p_vals=[]
    #Cycle through iterations of subsetting and calculating p-values
    for i in range(sample_size): 
        #Subset data
        ordered_arrays = ordered_subset(adata)
        random_arrays = random_subset(adata)
        #Calculate p-values for each gene and add as list to the p_val lists
        ordered_p_vals.append(execute_wilcox(ordered_arrays[0], ordered_arrays[1]))
        random_p_vals.append(execute_wilcox(random_arrays[0], random_arrays[1]))
    	#For every list of p-vals in the list of lists calculate number of signficant results given fdr correct at alpha input
        ordered_p_vals[-1].append(fdr_correct_sig(ordered_p_vals[-1], alpha))
        random_p_vals[-1].append(fdr_correct_sig(random_p_vals[-1], alpha))
        #Print update to console
        print("Iteration " + str(i) + " complete.")
    #write files
    with open(ordered_file, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in ordered_p_vals:
            writer.writerow(i)
        writeFile.close()
    with open(random_file, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in random_p_vals:
            writer.writerow(i)
        writeFile.close()

########################################################Subset Data Functions######################################################################################

'''
This function creates two approximately equal subsets of the original data in adata based on the order in the original AnnData object.  It assumes that the AnnData object is preordered by batch.  It ensures
even cluster representation assuming that biological differences are represented by cluster.  
'''
def ordered_subset(adata):
    
    #Get cluster and batch names
    cluster_names=list(set(adata.obs["cluster assignment"].tolist()))
    batch_names=list(set(adata.obs["batch"].tolist()))
    
    #Create objects to hold blank arrays for each of the two batches
    array_1 = adata.chunk_X(0)
    array_2 = adata.chunk_X(0)
    
    #Given the way we batch corrected the AnnData object is currently 
    #Cycle through clusters and divide up each cluster
    for i in cluster_names:
        
        #Store anndata object with just the cluster values in temp and get number of cells
        temp = adata[adata.obs["cluster assignment"] == i, :]
        temp.obs["value"] = 0
        cluster_n = temp.shape[0]
        
        #Figure out start position by random selection of values less than cluster_n
        start_pos=random.sample(range(cluster_n), 1)[0]
        #Get ordered indices
        subset_indices = get_ordered_indices(start_pos, cluster_n)
        #Pull first cluster_n//2 entries into one matrix and the remainig into another
        array_1 = np.vstack([array_1, temp.chunk_X(subset_indices)])
        array_2 = np.vstack([array_2, temp.chunk_X([x for x in range(cluster_n) if x not in subset_indices])])
    
    #Return both arrays
    return array_1, array_2


def get_ordered_indices(start_pos, total_values):
    
    #Calculate half total values
    half_total=total_values//2
    #Determine whether half_total + start_pos > total_values and proceed accordingly
    if start_pos + half_total <= total_values:
        #No need to circle back to beginning:
        out_indices = [x for x in range(half_total + start_pos) if x not in range(start_pos)]
        return out_indices
    else:
        #Here we have to circle back to pick up the rest of the indices
        out_indices = [x for x in range(total_values) if x not in range(start_pos)]
        remain_indices = [x for x in range(total_values) if x not in range(start_pos + half_total - total_values)]
        out_indices.extend(remain_indices)
        return out_indices

'''
This function creates two approximately equal subsets of the original data in adata based on random selection within each cluster.  It ensures even cluster representation assuming that biological differences are represented by cluster.  
'''
def random_subset(adata):
    
    #Get cluster and batch names
    cluster_names=list(set(adata.obs["cluster assignment"].tolist()))
    batch_names=list(set(adata.obs["batch"].tolist()))
    
    #Create objects to hold blank arrays for each of the two batches
    array_1 = adata.chunk_X(0)
    array_2 = adata.chunk_X(0)
    
    #Given the way we batch corrected the AnnData object is currently 
    #Cycle through clusters and divide up each cluster
    for i in cluster_names:
        
        #Store anndata object with just the cluster values in temp and get number of cells
        temp = adata[adata.obs["cluster assignment"] == i, :]
        temp.obs["value"] = 0
        cluster_n = temp.shape[0]
        
        #Get ordered indices
        subset_indices = random.sample(range(cluster_n),cluster_n//2)
        #Pull first cluster_n//2 entries into one matrix and the remainig into another
        array_1 = np.vstack([array_1, temp.chunk_X(subset_indices)])
        array_2 = np.vstack([array_2, temp.chunk_X([x for x in range(cluster_n) if x not in subset_indices])])
        
    #Return both arrays
    return array_1, array_2

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
        p_values.append(sp.stats.ranksums(array_1.T[i], array_2.T[i])[1])
    #Return p_values
    return p_values

def fdr_correct_sig(p_vals, alpha):
    #Fdr correct the list and return number below threshold
    fdrs = multi.fdrcorrection(p_vals, alpha)
    return sum(fdrs[0])

#########################################################Read in Arguments#############################################

if len(sys.argv) == 4:         
    #Call main thread
    main_thread(naming_conv=sys.argv[1], output_dir=sys.argv[2], input_file=sys.argv[3], sample_size=100, alpha=0.05)
elif len(sys.argv) == 6:
    #Call main thread function
    main_thread(naming_conv=sys.argv[1], output_dir=sys.argv[2], input_file=sys.argv[3], sample_size=sys.argv[4], alpha=sys.argv[5])
else:
    print("Incorrect number of parameters given.  Please try again.")
    print("Paramters should be the following: naming convention, output directory, input file(full path), sample size (optional/default=100), alpha (optional/default=0.05)")

