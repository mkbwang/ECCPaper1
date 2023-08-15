#!/bin/bash

#set input directory
in=/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/WGS/raw_data/
 
# get count of files in this directory
NUMFILES=$(ls -1 ${in}*_1.fastq.gz | wc -l)

# subtract 1 as we have to use zero-based indexing (first element is 0)
ZBNUMFILES=$(($NUMFILES - 1))

# submit array of jobs to SLURM
if [ $ZBNUMFILES -ge 0 ]; then
  sbatch --array=0-$ZBNUMFILES array_trim.pbs
else
  echo "No jobs to submit, since no input files in this directory."
fi
