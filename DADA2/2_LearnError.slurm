#!/bin/sh
#SBATCH --job-name=ECC1_LEARN_ERROR
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/learn_error_output.log
#SBATCH --error=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/learn_error_error.log


#load R

module load Rtidyverse/4.2.0

Rscript LearnError.R
