#!/bin/sh
#SBATCH --job-name=ECC1_16S_filter_trim
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --array=1-89
#SBATCH --mem-per-cpu=10g
#SBATCH --time=01:00:00
#SBATCH --output=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/filter_trim.log
#SBATCH --error=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/filter_trim_error.log

#load R

module load Rtidyverse/4.2.0

Rscript FilterAndTrim.R
