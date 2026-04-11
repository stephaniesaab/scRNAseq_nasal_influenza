#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8 #CPUS for parellization
#SBATCH --mem=128G
#SBATCH --job-name=singleR
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module load r-bundle-bioconductor/3.21

#Run my Rscript
Rscript 5.annotating_clusters.R
