#!/usr/bin/env Rscript

library(data.table, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)


f <-
m <- "/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/GCF_000001405.25_GRCh37.p13_rna.gbff.mRNA_ID-protein_ID.csv"
s <- "/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/Mutations.mutOut"
l <- "/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/GCF_000001405.25_GRCh37.p13_protein_ID-lengths.csv"

exonic <- fread(f, header=F, data.table=F)
map <- fread(m, header=F, data.table=F)
snap_exist <- fread(s, header=F, data.table=F)
protein_length <- fread(l, header=T, data.table=F)

colnames(map) <- c("Transcript", "ProteinID")
colnames(snap_exist) <- c("Transcript", "aaMutation", "SNAPscore")

# Add information to exonic file:
exonic_flted <- exonic %>% 
  filter(V3 != "UNKNOWN") %>%
  select(V1, V2, V3, V9) %>%
  separate_rows(V3, sep=",") %>%  # split one row into multiple by V3
  filter(nchar(V3)!=0) %>%
  separate(V3, c("Gene", "Transcript", "exon", "nuChange", "aaChange"), sep=":")  %>%
  mutate(aaMutation = gsub("p.", "", aaChange)) %>%
  left_join(map, by="Transcript") %>%
  left_join(protein_length, by="ProteinID") %>%
  left_join(snap_exist, by=c("Transcript", "aaMutation"))

# There are some transcripts in the ANNOVAR reference file that were deleted by RefSeq db
# So, there will be some (a small part) records unmapped
# unique(exonic_flted$Transcript[is.na(exonic_flted$ProteinID)])

exonic_mapped <- exonic_flted[!is.na(exonic_flted$ProteinID),]

exonic_mapped_unique <- exonic_mapped %>%
  group_by(V1) %>% top_n(1, Length)

any(duplicated(exonic_mapped_unique$V1))

# Originally, exonic does not have duplicated rows
# But after adding more information, there are duplicated rows in exonic_flted, because one variant maps to multiple transcript/proteins
# We take only the longest protein:
exonic_flted$ProtLen <- nchar(exonic_flted$sequence)
exonic_flted_nodup <- exonic_flted %>%
  group_by(V1) %>% 
  top_n(1, wt=ProtLen)

# Note that exonic_flted_nodup still has duplicated rows because isoforms are the same length

# Regardless, make some query files:
# (1) by uniprot ID:
query_1 <- exonic_flted_nodup[, c("UniProt.entry.ID", "aaChange")]
query_1$aaChange <- gsub("p.", "", query_1$aaChange)
query_1 <- query_1[!is.na(query_1$UniProt.entry.ID),]
write.table(query_1, "/Users/WangYanran/Desktop/Rerun_Crohns/Step3/query_UniProt-ID.txt",
            quote=F, sep="\t", row.names=F, col.names=F)
# (2) by sequence:
query_2 <- exonic_flted_nodup[, c("UniProt.entry.ID", "aaChange")]


# %>%
#   filter(!is.na(Prot_length)) %>%
#   mutate(Mut=gsub("p.", "", aaChange)) %>%
#   left_join(snap_score, by=c("Transcript", "Mut"))