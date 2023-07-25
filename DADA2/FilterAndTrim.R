library(dada2); packageVersion("dada2")


# File parsing 
pathF <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/FWD" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/REV" # CHANGE ME ...
filtpathF <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/FWD" # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/REV" # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
subset <- seq((ID-1)*10+1, min(ID*10, length(fastqFs)))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
RunOut<-filterAndTrim(fwd=file.path(pathF, fastqFs[subset]), filt=file.path(filtpathF, fastqFs[subset]),
              rev=file.path(pathR, fastqRs[subset]), filt.rev=file.path(filtpathR, fastqRs[subset]),
              trimLeft = c(9, 9), truncLen=c(240,210), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

save(RunOut, file="/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/filter_trim_dada2.Rdata")

