#!/bin/bash
#SBATCH --job-name=mcore_job
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --ntasks=1

module load R/4.3.2

Rscript --no-save R/Chapter_02_standalone.R
