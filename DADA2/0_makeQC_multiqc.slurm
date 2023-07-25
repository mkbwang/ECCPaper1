#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=make_multi_qc_ECCPaper1
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=BEGIN,END,ERROR
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB 
#SBATCH --time=10:00:00
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/qc_output.log
#SBATCH --error=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/DADA2/qc_error.log

module load python3.9-anaconda
conda activate microbiome # Mukai's local environment
module load Bioinformatics
module load fastqc

## forward
cd /nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/FWD
fastqc *.fastq.gz --quiet --outdir=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/FWD/fastqcresults
multiqc ./ --quiet --force --filename  FWD_16S_multiqc_report.html --outdir /nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/qc_summary/


## backward
cd /nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/REV
fastqc *.fastq.gz --quiet --outdir=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/REV/fastqcresults
multiqc ./ --quiet --force --filename  REV_16S_multiqc_report.html --outdir /nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/qc_summary/



