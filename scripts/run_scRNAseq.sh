#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8 #8 CPUS for parellization
#SBATCH --mem=64G
#SBATCH --job-name=scRA_process
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module purge
module load r-bundle-bioconductor/3.21

#Run my Rscript
