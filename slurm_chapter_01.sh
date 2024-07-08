#!/bin/bash
#SBATCH --job-name=run_method_experiments
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=32GB
module load R/4.3.2

Rscript --no-save R/chapter_01_samplingonly.R
