library(dada2); packageVersion("dada2")
library(ggplot2)

# File parsing 
filtpathF <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/FWD" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/REV" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), '_'),`[`,1) # Assumes filename = samplename_XXX.fastq.gz # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#Plot & save errors 
Run_errF<-plotErrors(errF, nominalQ=TRUE)
Run_errR<-plotErrors(errR, nominalQ=TRUE)
ggsave(filename=paste(filtpathF, "run_errorFwd.pdf", sep="/"), plot=Run_errF)
ggsave(filename=paste(filtpathR, "run_errorRev.pdf", sep="/"), plot=Run_errR)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
dadaFs<-vector("list", length(sample.names))
dadaRs<-vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaFs)<-sample.names
names(dadaRs)<-sample.names
for(sam in sample.names) {
  # cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# Construct sequence table and remove chimeras
seqtab_run <- makeSequenceTable(mergers)
saveRDS(seqtab_run, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/seqtab_16S.rds")# CHANGE ME to where you want sequence table saved
#Sanity check -- loss through steps
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)
track_run1<-track
rm(derepF); rm(derepR)
