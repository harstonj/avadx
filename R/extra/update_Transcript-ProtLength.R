#!/usr/bin/env Rscript
library(seqinr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

dbfolder <- args[1]

setwd(dbfolder)
#df <- read.csv("Transcript-ProtLength.csv", stringsAsFactors=F, header=T)

fas <- read.fasta("prot_seqs.txt", as.string=T, seqtype="AA", set.attributes=F, whole.header=T)
fas_names <- names(fas)
fas_mRNA <- strsplit(fas_names, "\\.|\\|")
fas_mRNA <- unlist(lapply(fas_mRNA, function(x) x[2]))

gene_names <- regmatches(fas_names, gregexpr("\\[.+?\\]", fas_names))
gene_names <- lapply(gene_names, function(x) substr(x, 2, nchar(x)-1))
gene_names <- lapply(gene_names, function(x) x[startsWith(x, "gene=")])
gene_names <- lapply(gene_names, function(x) gsub("gene=", "", x))
gene_names <- unlist(gene_names)

fas_lengths <- c()
for(i in 1:length(fas)){
  fas_lengths <- c(fas_lengths, nchar(fas[[i]]))
}

df <- data.frame("Gene"=gene_names,
                 "Transcript"=fas_mRNA,
                 "Prot_length"=fas_lengths)

setwd(dbfolder)
write.csv(df, "Transcript-ProtLength_update.csv", quote=F, row.names=F)
