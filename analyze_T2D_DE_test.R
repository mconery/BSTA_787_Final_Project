#This script analyzes the results of the T2D DE tests.  It exports a histogram comparing the number of DE genes for 
#each iteration of the ordered testing method.  Assumes the output files from the DE tests are al stored in the input directory.  


# 0) Load libraries ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(metap)
library(car)

# 1) Read in data and Minimally Process ====

#Set directories and alpha
alpha <- 0.05
#Set Input Directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/processed_data/"
#Set output directory
out_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/processed_data/"


#########################Cell Types##########################

#Load the datasets
raw_results <- read.table(paste0(inp_dir, "T2D_DE_cell_type_results.tsv"), sep = "\t", header = TRUE, row.names = 1)

#Create vectors to be cluster names and number of sig results
clusters <- vector()
sig_results <- vector()
#Make copy for fdr corrections
fdr_results <- raw_results
#Cycle through clusters and fdr correct p_vals
for (i in 1:(ncol(fdr_results)/2)) {
  #overwrite column data
  fdr_results[,2*i] <- p.adjust(fdr_results[,2*i], method = "fdr")
  #Append info to vectors
  clusters <- append(clusters, paste0("batch_", substr(colnames(fdr_results)[2*i], 1, gregexpr(pattern ='_', colnames(fdr_results)[2*i])[[1]][1] - 1)))
  sig_results <- append(sig_results, sum(fdr_results[,2*i] < alpha))
}

#Cbind vectors and export table
export_table <- cbind.data.frame(cluster=clusters, DE_genes=sig_results, Pct_DE=sig_results/nrow(fdr_results))
write.table(export_table, file = paste0(out_dir, "T2D_DE_cell_type_Genes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

##########################Clusters##################################

#Load the datasets
raw_results <- read.table(paste0(inp_dir, "T2D_DE_cluster_results.tsv"), sep = "\t", header = TRUE, row.names = 1)

#Create vectors to be cluster names and number of sig results
clusters <- vector()
sig_results <- vector()
#Make copy for fdr corrections
fdr_results <- raw_results
#Cycle through clusters and fdr correct p_vals
for (i in 1:(ncol(fdr_results)/2)) {
  #overwrite column data
  fdr_results[,2*i] <- p.adjust(fdr_results[,2*i], method = "fdr")
  #Append info to vectors
  clusters <- append(clusters, paste0("batch_", substr(colnames(fdr_results)[2*i], 1, gregexpr(pattern ='_', colnames(fdr_results)[2*i])[[1]][1] - 1)))
  sig_results <- append(sig_results, sum(fdr_results[,2*i] < alpha))
}

#Cbind vectors and export table
export_table <- cbind.data.frame(cluster=clusters, DE_genes=sig_results, Pct_DE=sig_results/nrow(fdr_results))
write.table(export_table, file = paste0(out_dir, "T2D_DE_cluster_Genes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
