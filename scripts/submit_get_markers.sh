#!/bin/bash
#SBATCH --job-name=Marker_stable
#SBATCH --account=def-itobias
#SBATCH --time=02:00:00             
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1         
#SBATCH --mem=128G                  # for 154k cells
#SBATCH --output=%x-%j.out

# Load the R module (check which versions are available on Narval)
module load r-bundle-bioconductor/3.21

# Run the script
Rscript 4.checking_cluster_csv.R
