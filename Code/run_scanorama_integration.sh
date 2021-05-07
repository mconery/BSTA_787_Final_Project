#!/bin/bash

#SBATCH -o run_scanorama_integration.log
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

#Activate environment
source activate scanoramaenv

#Execute python program
python ~/BSTA_787/scanorama_correct.py
