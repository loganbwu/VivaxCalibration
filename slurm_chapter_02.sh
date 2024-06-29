#!/bin/bash
#SBATCH --job-name=mcore_job
#SBATCH --partition=long
#SBATCH --time=4-00
#SBATCH --cpus-per-task=96
#SBATCH --mem=120G
#SBATCH --ntasks=1
#SBATCH --mail-user=wu.l@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R/4.4

Rscript --no-save R/Chapter_02_standalone.R
