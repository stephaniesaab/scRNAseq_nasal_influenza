#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1 #Only need one CPU because running sequential
#SBATCH --mem=192G #High RAM to handle residual matrix of SCTransform
#SBATCH --job-name=umap_processing
#SBATCH --output=%x-%j.out #Log file for errors or output

#load modules
module load r-bundle-bioconductor/3.21

#Run my Rscript
Rscript 2.cluster_annotation_plot.R
