#!/bin/bash
#SBATCH --job-name=chapter_02_long
#SBATCH --partition=long
#SBATCH --time=14-00
#SBATCH --cpus-per-task=96
#SBATCH --mem=256G
#SBATCH --ntasks=1
#SBATCH --mail-user=wu.l@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

max_hours=240

module load R/4.4

Rscript --no-save R/Chapter_02_standalone.R "$max_hours"
