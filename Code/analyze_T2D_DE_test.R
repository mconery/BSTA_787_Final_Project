#This script analyzes the results of the T2D DE tests.  It exports lists showing the number of DE genes for 
#each cell type and cluster.  It also makes a scatterplot of number of significant DE vs cluster size.
#Assumes the output files from the DE tests are all stored in the input directory.  


# 0) Load libraries ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(metap)
library(car)

# 1) Read in data and Make lists ====

#Set directories and alpha
alpha <- 0.05
#Set Input Directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/processed_data/"
#Set output directories
out_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/processed_data/"
plot_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/figures/"


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
  clusters <- append(clusters, paste0(substr(colnames(fdr_results)[2*i], 1, gregexpr(pattern ='_', colnames(fdr_results)[2*i])[[1]][1] - 1)))
  sig_results <- append(sig_results, sum(fdr_results[,2*i] < alpha))
}

#Cbind vectors and export table
export_table <- cbind.data.frame(cell_type=clusters, DE_genes=sig_results, Pct_DE=sig_results/nrow(fdr_results))
write.table(export_table, file = paste0(out_dir, "T2D_DE_cell_type_Genes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Store cell types
cell_type_df <- export_table[,c("cell_type", "DE_genes")]
colnames(cell_type_df)[1] <- "type"

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
  clusters <- append(clusters, paste0("cluster_", substr(colnames(fdr_results)[2*i], 2, gregexpr(pattern ='_', colnames(fdr_results)[2*i])[[1]][1] - 1)))
  sig_results <- append(sig_results, sum(fdr_results[,2*i] < alpha))
}

#Cbind vectors and export table
export_table <- cbind.data.frame(cluster=clusters, DE_genes=sig_results, Pct_DE=sig_results/nrow(fdr_results))
write.table(export_table, file = paste0(out_dir, "T2D_DE_cluster_Genes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 2) Create scatterplot of cluster/cell type size and # of sig genes====

#Read in obs
obs <- read.table(paste0(inp_dir, "obs.csv"), sep = ",", header = TRUE, row.names = 1)
#Remove nans
obs <- obs[which(obs$cell.type != ''),]
#Sort observations
obs <- obs[order(obs$cell.type),]
#Make table of observations
obs_table <- table(obs$cell.type, obs$cluster.assignment)

#Remove nan from cell type df and sort df
cell_type_df <- cell_type_df[which(cell_type_df$type != "nan"),]
cell_type_df <- cell_type_df[order(cell_type_df$type),]
#Bind on number of cells and grouping
cell_type_df <- cbind.data.frame(cell_type_df, num_cells=rowSums(obs_table), grouping=rep("cell type", nrow(cell_type_df)))
#Calculate spearman correlation
cell_type_cor <- cor(cell_type_df$DE_genes, cell_type_df$num_cells, method = "spearman")

#Save export_table as cluster df
cluster_df <- export_table
#Remove "cluster" 
cluster_df$cluster <- as.numeric(substr(cluster_df$cluster, 9, nchar(cluster_df$cluster)))
#Remove any clusters that don't get tested and sort
cluster_df <- cluster_df[which(as.character(cluster_df$cluster) %in% colnames(obs_table)),]
cluster_df <- cluster_df[order(cluster_df$cluster),]
colnames(cluster_df)[1] <- "type" 
#Append column to cluster_df
cluster_df <- cbind.data.frame(cluster_df, num_cells=colSums(obs_table))
#Drop pct column and add column for grouping
cluster_df <- cluster_df[,c(1,2,4)]
cluster_df <- cbind.data.frame(cluster_df, grouping=rep("cluster", nrow(cluster_df)))
#Calculate correlation
cluster_cor <- cor(cluster_df$DE_genes, cluster_df$num_cells, method = "spearman")

#Make master_df
master_df <- rbind.data.frame(cell_type_df, cluster_df)

#Make Scatterplots
tiff(paste0(plot_dir, "correlation_scatterplot.jpg"), width = 11.5, height = 8, units = 'in', res = 500)
ggplot(master_df, aes(x=num_cells, y=DE_genes, color=grouping, shape=grouping)) +
  geom_point() + geom_smooth(method=lm, se = FALSE) + 
  xlab("Number of Cells per Group") + ylab("DE Genes Detected") + 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  theme(panel.background = element_blank()) + 
  annotate("text", x=700, y=12100, label = paste0("cluster spearman correlation = ", format(cluster_cor, digits=3)), color="#00BFC4", hjust="0", size = 5) + 
  annotate("text", x=3000, y=0, label = paste0("cell type spearman correlation = ", format(cell_type_cor, digits=3)), color="#F8766D", hjust="0", size = 5)
dev.off()
  




