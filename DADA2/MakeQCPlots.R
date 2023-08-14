library(dada2); packageVersion("dada2")
library(ggplot2)
####
resultpath<- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/qcplots"
################RUN 1 ###################################################################################################
pathF <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/FWD/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/raw_data/REV/" # CHANGE ME ...
fastqFs <- paste(pathF, sort(list.files(pathF, pattern="*.fastq.gz")), sep="")
fastqRs <- paste(pathR, sort(list.files(pathR, pattern="*.fastq.gz")), sep="")
#plot forwards
max <- 9
x <- seq_along(fastqFs)
FwdList<-split(fastqFs, ceiling(x/max))
FWDplots <- lapply(FwdList, plotQualityProfile)
names(FWDplots)=seq(1:length(FwdList))
lapply(names(FWDplots), 
       function(x) ggsave(filename=paste(resultpath, "Forward_", x, ".pdf", sep=""), device="pdf", plot=FWDplots[[x]]))
#plot reverse
x <- seq_along(fastqRs)
RevList<-split(fastqRs, ceiling(x/max))
Revplots <- lapply(RevList, plotQualityProfile)
names(Revplots)=seq(1:length(RevList))
lapply(names(Revplots), 
       function(x) ggsave(filename=paste(resultpath, "Reverse_", x, ".pdf", sep=""), device="pdf", plot=Revplots[[x]]))

