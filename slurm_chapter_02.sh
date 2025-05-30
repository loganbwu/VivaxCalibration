#!/bin/bash
#SBATCH --job-name=chap_02
#SBATCH --partition=regular
#SBATCH --time=2-00
#SBATCH --cpus-per-task=128
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --mail-user=wu.l@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

max_hours=24

module load R/4.4

Rscript --no-save R/Chapter_02_standalone.R "$max_hours"
