#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4 #4 CPUS for parellization
#SBATCH --mem=128G
#SBATCH --job-name=DE_umap
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module load r-bundle-bioconductor/3.21

#Run my Rscript
Rscript 3.DE_umap.R
