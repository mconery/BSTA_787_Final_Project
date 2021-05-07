#!/bin/bash

#SBATCH -o run_DE_scvi.log
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

#Activate environment
source activate scvienv

#Define needed locations
output_dir=~/BSTA_787/processed_pancreas_data/scvi
input_file=$output_dir/scvi_corrected.h5ad

#Execute python program
python ~/BSTA_787/test_DE.py scvi $output_dir $input_file
