rm(list=ls())
library(dplyr)
folder <- "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1"

# run on command line:  esearch -db sra -query PRJNA752888 | efetch -format runinfo > runinfo_all.csv
SRR_metadata <- read.csv(file.path(folder, "metadata", "runinfo_all.csv"))


WGS_metadata <- SRR_metadata %>% filter(LibraryStrategy == "WGS")
Amplicon16S_metadata <- SRR_metadata %>% filter(LibraryStrategy == "AMPLICON")

write.csv(WGS_metadata, file.path(folder, 'WGS', 'runinfo_WGS.csv'),
          row.names=FALSE, quote=FALSE)

write.csv(Amplicon16S_metadata, file.path(folder, 'Amplicon16S', 'runinfo_Amplicon16S.csv'),
          row.names=FALSE, quote=FALSE)

WGSmeta_subset <- WGS_metadata %>% select(Run, SampleName)
Amplicon16S_subset <- Amplicon16S_metadata %>% select(Run, SampleName) %>%
  rename(FoxCavID=SampleName)




