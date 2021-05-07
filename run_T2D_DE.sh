#!/bin/bash

#SBATCH -o run_T2D_DE.log
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

#Activate environment
source activate scanoramaenv

#Define needed locations
output_dir=~/BSTA_787/processed_pancreas_data/scanorama
input_file=$output_dir/scanorama_corrected.h5ad
T2D_annotations=/home/conerym/BSTA_787/raw_pancreas_data/diabetes_annotations.tsv

#Execute python program
python ~/BSTA_787/DE_T2D.py $output_dir $input_file $T2D_annotations


