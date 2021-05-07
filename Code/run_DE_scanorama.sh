#!/bin/bash

#SBATCH -o run_DE_scanorama.log
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

#Activate environment
source activate scanoramaenv

#Define needed locations
output_dir=~/BSTA_787/processed_pancreas_data/scanorama
input_file=$output_dir/scanorama_corrected.h5ad

#Execute python program
python ~/BSTA_787/test_DE.py scanorama $output_dir $input_file


