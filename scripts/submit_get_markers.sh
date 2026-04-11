#!/bin/bash
#SBATCH --job-name=Marker_stable
#SBATCH --account=def-itobias
#SBATCH --time=02:00:00             # 2 hours is plenty with 16 CPUs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1         # Matches the workers in the R script
#SBATCH --mem=128G                  # Safe RAM for 154k cells
#SBATCH --output=%x-%j.out

# Load the R module (check which versions are available on Narval)
module load r-bundle-bioconductor/3.21

# Run the script
Rscript get_markers0.4.R
