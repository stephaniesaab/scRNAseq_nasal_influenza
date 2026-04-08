#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2 #2 CPUS for parellization
#SBATCH --mem=32G
#SBATCH --job-name=filter_cells
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module load r-bundle-bioconductor/3.21

#Run my Rscript
Rscript 1.filtering_scRNA.R
