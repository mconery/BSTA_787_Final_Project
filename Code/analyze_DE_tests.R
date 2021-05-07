#This script analyzes the results of the iterated DE tests.  It exports a histogram comparing the number of DE genes for 
#each iteration of the ordered testing method.  Assumes the output files from the DE tests are al stored in the input directory.  

# 0) Read in Needed Arguments from Command ====

#Set jump variable representing junk added to args (This has been determined by trial and error)
jump = 5

#Need to read in inp_dir and out_dir
args <- commandArgs();
if (length(args) == jump + 2) {
  
  #Print confirmation output
  print("Necessary Parameters Input")
  
  #Threshold for Expression (Because of the conversion to TPM the expression levels should be independent of the sequencing depth)
  inp_dir <- as.character(args[jump + 1])
  out_dir <- as.character(args[jump + 2])
  
} else {
  print("Usage: %> Rscript analyze_DE_tests.R inp_dir out_dir");
  quit(save="no");
}

#Testing locs
##Set Input Directory
#inp_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/processed_data/"
##Set output directory
#out_dir <- "C:/Users/mitch/Documents/UPenn/BSTA 787/Final_Project/figures/"

# 1) Load libraries ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(metap)
library(car)

# 2) Read in data and Minimally Process ====

#Load the datasets
scanorama_order <- read.table(paste0(inp_dir, "scanorama_ordered_pvals.tsv"), sep = "\t", header = FALSE)
scanorama_random <- read.table(paste0(inp_dir, "scanorama_random_pvals.tsv"), sep = "\t", header = FALSE)
scvi_order <- read.table(paste0(inp_dir, "scvi_ordered_pvals.tsv"), sep = "\t", header = FALSE)
scvi_random <- read.table(paste0(inp_dir, "scvi_random_pvals.tsv"), sep = "\t", header = FALSE)

#Get number of significant results per iteration
scanorama_order_sig <-  scanorama_order[,ncol(scanorama_order)]
scanorama_random_sig <-  scanorama_random[,ncol(scanorama_random)]
scvi_order_sig <-  scvi_order[,ncol(scvi_order)]
scvi_random_sig <-  scvi_random[,ncol(scvi_random)]

# 3) Process data and return desired histogram ====

#Make dataframe of sig_results
sig_results = cbind.data.frame(correction_method = append(rep("scanorama", length(scanorama_order_sig)), rep("scvi", length(scvi_order_sig))), 
                               num_sig=append(scanorama_order_sig, scvi_order_sig))
#Calculate categorical averages
med <- ddply(sig_results, "correction_method", summarise, grp.median=median(num_sig))
# Add mean lines
tiff(paste0(out_dir, "comparison_histogram.jpg"), width = 11.5, height = 8, units = 'in', res = 500)
ggplot(sig_results, aes(x=num_sig, fill=correction_method)) + 
  geom_histogram(position = "identity", alpha = 0.2, binwidth = 100) +
  geom_vline(data=med, aes(xintercept=grp.median, color=correction_method), linetype="dashed") +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  xlab("DE Genes per Iteration") + ylab("Number of Iterations")
dev.off()

