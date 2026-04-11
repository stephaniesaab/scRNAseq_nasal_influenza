#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8 #CPUS for parellization
#SBATCH --mem=128G
#SBATCH --job-name=DE_GSEA
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module load r-bundle-bioconductor/3.21

#Run my Rscript
Rscript 6.DE_GSEA.R
