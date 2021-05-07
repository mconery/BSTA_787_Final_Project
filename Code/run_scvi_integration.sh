#!/bin/bash

#SBATCH -o run_scvi_integration.log
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

#Activate environment
source activate scvienv

#Execute python program
python ~/BSTA_787/scvi_correct.py
