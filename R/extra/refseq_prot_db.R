library(seqinr, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)

p <- "./db/GCF_000001405.25_GRCh37.p13_protein.faa.gz"
proteins <- read.fasta(p, as.string=T, seqtype="AA")

proteins_len <- data.frame("ID"=rep(NA, length(proteins)), "Length"=rep(NA, length(proteins)))
for(i in 1:length(proteins)){
  print(i)
  proteins_len$ID[i] <- names(proteins)[i]
  proteins_len$Length[i] <- nchar(proteins[[i]][1])
}

proteins_len$ProteinID <- unlist(lapply(strsplit(proteins_len$ID, "\\."), function(x) x[1]))
proteins_len$ID <- NULL

# Write complete lengths info:
write.csv(proteins_len[,c(2,1)], "./db/GCF_000001405.25_GRCh37.p13_protein_ID-lengths.csv", quote=F, row.names=F)

# Remove if lengths are the same:
m <- "./parse_result_protein.csv"
gp <- data.table::fread(m, data.table=F, header=F)
colnames(gp) <- c("ProteinID", "Gene")

gp <- merge(gp, proteins_len, by="ProteinID", all=T)

gp_unique <- gp %>%
  group_by(Gene) %>%
  top_n(1, Length)

#sum(duplicated(gp_unique$Gene))
#View(gp_unique[duplicated(gp_unique$Gene),])

gp_unique$ProteinID.numeric <- as.numeric(substr(gp_unique$ProteinID, 4, nchar(gp_unique$ProteinID)))

gp_unique2 <- gp_unique %>%
  group_by(Gene) %>%
  top_n(-1, ProteinID.numeric) %>%
  select(-ProteinID.numeric)

sum(duplicated(gp_unique2$Gene))
# Write unique lengths:
write.csv(gp_unique2, "./db/GCF_000001405.25_GRCh37.p13_protein_ID-lengths_unique.csv", quote=F, row.names=F)

