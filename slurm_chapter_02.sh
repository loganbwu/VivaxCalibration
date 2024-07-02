#!/bin/bash
#SBATCH --job-name=mcore_job
#SBATCH --partition=regular
#SBATCH --time=4-00
#SBATCH --cpus-per-task=128
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --mail-user=wu.l@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R/4.4

Rscript --no-save R/Chapter_02_standalone.R
