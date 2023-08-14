library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st_16S <- readRDS("/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/seqtab_16S.rds")

# Remove chimeras
seqtab <- removeBimeraDenovo(st_16S, method="consensus", multithread=TRUE)
dim(st_16S)
dim(seqtab)
sum(st_16S)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- rowSums(seqtab)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track_postchimera<-track
saveRDS(track_postchimera, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/track_postchimera_HOMD.rds")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/HOMD_assigntaxa_togenus_plusextras.fasta", multithread=TRUE)
taxa <- addSpecies(taxa, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/HOMD_AssignSpecies_plusextras.fasta")

taxa_toSp <- assignTaxonomy(seqtab, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/HOMD_AssignTaxaToSpecies_plusextras.fasta", multithread = TRUE)

# Write to disk
saveRDS(seqtab, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/seqtab_final_homd.rds") # CHANGE ME to where you want sequence table saved
saveRDS(taxa, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/tax_final_homd.rds") # CHANGE ME ...
saveRDS(taxa_toSp, "/nfs/turbo/sph-bfoxman/People/wangmk/ECCPaper1/Amplicon16S/filtered/tax_toSp_final_homd.rds") # CHANGE ME ...

