#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --job-name=ECCdownload_WGS
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mem-per-cpu=1g
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/WGS/download_WGS.log
#SBATCH --error=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/WGS/download_WGS_error.log

folder=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1

module load Bioinformatics
module load gcc
module load sratoolkit

sed 1d runinfo_WGS.csv | cut -d "," -f 1 > SRR.numbers.wgs
# cat SRR.numbers.wgs | parallel prefetch 
cat SRR.numbers.wgs | parallel fastq-dump --gzip --readids --split-files --dumpbase --skip-technical -O ${folder}/WGS/raw_data
