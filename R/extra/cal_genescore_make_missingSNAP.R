#!/usr/bin/env Rscript

## This script is for gene_score calculation.
## Input and output of this script is only for one .exonic_variant_function file, i.e. one individual/sample
## For calculation of multiple samples, user can either:
##  (1) Use a for loop in the shell script which calls this Rscript.
##  (2) Use cluster (amarel) and run this paralelly for each input files using arrays.

library(optparse, quietly=T, warn.conflicts=F)
library(data.table, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)

# Example inputs:
# f <- "/path/to/ANNOVAR/output/files/CD_01.exonic"
# snapscore_file <- "/path/to/db/folder/Mutations.mutOut"
# out <- "/path/to/output/folder"

# f <- "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/exonic_variant_function/sample.25000.fa.avinput.exonic_variant_function"
# snapscore_file <- "/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/Mutations.mutOut"
# protlength_file <- "/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/Transcript-ProtLength_cleaned.csv"


option_list = list(
  make_option(c("-f", "--exonic_file"), type="character", default=NULL, 
              help=".exonic_variant_function by ANNOVAR of one individual", metavar="character"),
  make_option(c("-s", "--snapscore_file"), type="character", default=NULL, 
              help="snapscore file where columns are: gene (NM_ number), aa_mutation, snapscore", metavar="character"),
  make_option(c("-l", "--protlength_file"), type="character", default=NULL, 
              help="protein length file where columns are: gene name, NM_ number, corresponding protein length from RefSeq", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="path to the output folder", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

# Read-in SNAP scores:
snap_score <- fread(opt$snapscore_file, data.table=F)
colnames(snap_score) <- c("Transcript", "Mut", "SNAP_Score")

# Read-in Acc-protein_length:
protL <- fread(opt$protlength_file, data.table=F)
protL$Gene <- NULL

# Pre-precessing of the input file:
exonic <- fread(opt$exonic_file, data.table=F)
exonic_flted <- exonic %>% 
  filter(V3 != "UNKNOWN") %>%
  select(V1, V2, V3, V9) %>%
  separate_rows(V3, sep=",") %>%
  filter(nchar(V3)!=0) %>%
  separate(V3, c("Gene", "Transcript", "exon", "nuChange", "aaChange"), sep=":") %>%
  left_join(protL, by="Transcript") %>%
  filter(!is.na(Prot_length)) %>%
  mutate(Mut=gsub("p.", "", aaChange)) %>%
  left_join(snap_score, by=c("Transcript", "Mut"))

# Output missing SNAPs:
exonic_missingSNAP <- exonic_flted %>% 
  filter(V2=="nonsynonymous SNV" & is.na(SNAP_Score))
setwd(opt$out)
dir.create("missingSNAP")
for(transc in unique(exonic_missingSNAP$Transcript)){
  tmp <- exonic_missingSNAP[exonic_missingSNAP$Transcript==transc,]
  write.table(tmp$Mut, paste0("missingSNAP/", transc, ".mutation"), quote=F, col.names=F, row.names=F)
  rm(tmp)
}
