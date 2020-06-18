#!/usr/bin/env Rscript
library(dplyr)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)

exonic.file <- args[1]
dbfolder <- args[2]
outputfolder <- args[3]

setwd(dbfolder)
df_len <- read.csv("Transcript-ProtLength.csv", stringsAsFactors=F, header=T)
#mutOut <- read.table("Mutations.mutOut", stringsAsFactors=F, header=F, sep="\t")

exonic <- data.table::fread(exonic.file, header=F, data.table=F)
exonic <- exonic[exonic$V3!="UNKNOWN",]
exonic <- exonic[exonic$V2 %in% c("synonymous SNV", "nonsynonymous SNV", "stopgain", "stoploss"), ]

variants <- strsplit(exonic$V3, ",")
variants <- lapply(variants, function(x) strsplit(x, ":"))

variant_lengths <- lapply(variants, function(x) unlist(lapply(x, length)))
if(unique(unlist(variant_lengths))!=5){
  stop(".exonic_variant_function file has rows with weird rows, please check and run this script again.")
}

variants_mRNA <- lapply(variants, function(x) lapply(x, function(y) y[2]))
variants_aaMut <- lapply(variants, function(x) lapply(x, function(y) y[5]))
variants_gene <- lapply(variants, function(x) lapply(x, function(y) y[1]))

gene_mRNA_list <- list()
for(i in 1:length(variants)){
  g <- unlist(variants_gene[[i]])
  line <- exonic$V1[i]
  trans <- unlist(variants_mRNA[[i]])
  aaMut <- unlist(variants_aaMut[[i]])
  dfi <- data.frame("Gene"=g, "Transcript"=trans, "aaMut"=gsub("p.", "", aaMut), "line"=line)
  gene_mRNA_list[[i]] <- dfi
  rm(g, trans, dfi)
}

gene_mRNA <- data.table::rbindlist(gene_mRNA_list, use.names=T)
gene_mRNA$prot_length <- df_len$Prot_length[match(gene_mRNA$Transcript, df_len$Transcript)]

setwd(outputfolder)
  write.table(unique(as.character(gene_mRNA$Transcript[is.na(gene_mRNA$prot_length)])), 
              file="missing_mRNA_accession_global.txt", quote=F, col.names=F, row.names=F)
if(any(is.na(gene_mRNA$prot_length))){
  print(
    paste0(
      "Transcript-ProtLength.csv lacks protein length info for ",
      length(which(is.na(gene_mRNA$prot_length))), " GLOBAL transcripts (",
      round(length(which(is.na(gene_mRNA$prot_length)))/nrow(gene_mRNA)*100), "%). "
    )
  )
}

# check if Mutations.mutOut already contains the variants:
#gene_mRNA$Exist <- ifelse(paste0(gene_mRNA$Transcript, "-", gene_mRNA$aaMut) %in% paste0(mutOut$V1, "-", mutOut$V2), T, F)
#print(paste0(sum(gene_mRNA$Exist), " mutations (", round(sum(gene_mRNA$Exist)/nrow(gene_mRNA)*100),"%) already have SNAP scores in Mutations.mutOut. Another ", sum(!gene_mRNA$Exist), " new mutations need to run with SNAP."))
#print("Generating SNAP input files...")

# For the same line, keep only the mutation in the longest protein:
#gene_mRNA <- gene_mRNA[gene_mRNA$Exist==F, ]
gene_mRNA_nodup <- gene_mRNA %>%
  group_by(line) %>%
  filter(prot_length == max(prot_length))

setwd(outputfolder)
write.table(unique(as.character(gene_mRNA_nodup$Transcript[is.na(gene_mRNA_nodup$prot_length)])), 
            file="missing_mRNA_accession_query.txt", quote=F, col.names=F, row.names=F)
if(any(is.na(gene_mRNA_nodup$prot_length))){
  print(
    paste0(
      "Transcript-ProtLength.csv lacks protein length info for ",
      length(which(is.na(gene_mRNA_nodup$prot_length))), " QUERY transcripts (",
      round(length(which(is.na(gene_mRNA_nodup$prot_length)))/nrow(gene_mRNA_nodup)*100), "%)"
    )
  )
}

# Read in fasta files:
setwd(dbfolder)
fas <- read.fasta("prot_seqs.txt", as.string=T, seqtype="AA", set.attributes=F)
fas_names <- names(fas)
fas_mRNA <- strsplit(fas_names, "\\.|\\|")
fas_mRNA <- unlist(lapply(fas_mRNA, function(x) x[2]))

# 1. check if the prot_seq.txt contains the protein sequence:
missing_prot_seq <- setdiff(unique(as.character(gene_mRNA$Transcript)), fas_mRNA)

if(length(missing_prot_seq)>0){
  setwd(outputfolder)
  write.table(missing_prot_seq, "TranscriptAccess_missing_prot_seq_global.txt", quote=F, col.names=F, row.names=F)
  print(paste0(length(missing_prot_seq), " GLOBAL transcripts do not have protein sequences in the db folder. Please append the protein sequence into the prot_seqs.txt file! mRNA transcript accession numbers were output to the output folder and is named: TranscriptAccess_missing_prot_seq.txt. Please use NCBI batch entrez query to obtain protein sequences (https://www.ncbi.nlm.nih.gov/sites/batchentrez)."))
}

# # 2. make SNAP input:
# setwd(outputfolder)
# for(i in 1:nrow(gene_mRNA_nodup)){
#   aaMutation <- as.character(gene_mRNA_nodup$aaMut)[i]
#   pos <- as.numeric(substr(aaMutation, 2, nchar(aaMutation)-1))
#   ref <- substr(aaMutation, 1, 1)
#   alt <- substr(aaMutation, nchar(aaMutation), nchar(aaMutation))
#   trans <- as.character(gene_mRNA_nodup$Transcript)[i]
  
#   aaSeq <- fas[[match(trans, fas_mRNA)]]
#   # check if REF is in line with aaSeq: otherwise ignore this mutation
#   if(ref==substr(aaSeq, pos, pos) & alt!="X"){
#     write.table(aaMutation, paste0(trans, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
#   }
#   if(!file.exists(paste0(trans, ".fasta"))){
#     write.fasta(aaSeq, fas_names[match(trans, fas_mRNA)], file.out=paste0(trans, ".fasta"))
#   }
# }

# 1. check if prot_seq.txt contains all sequence for relevant query subset:
missing_prot_seq_nodup <- setdiff(unique(as.character(gene_mRNA_nodup$Transcript)), fas_mRNA)
setwd(outputfolder)
  write.table(missing_prot_seq_nodup, "TranscriptAccess_missing_prot_seq_query.txt", quote=F, col.names=F, row.names=F)

if(length(missing_prot_seq_nodup)>0){
  print(paste0(length(missing_prot_seq_nodup), " QUERY transcripts do not have protein sequences in the db folder"))
}

# 2. Generate snapcachedb input
setwd(outputfolder)
matched = match(unique(gene_mRNA_nodup$Transcript), fas_mRNA)
write.fasta(fas[matched], fas_mRNA[matched], file.out='snapfun_query.fa')
write.table(gene_mRNA_nodup[c('Transcript', 'aaMut')], "snapfun_query.ids", quote=F, col.names=F, row.names=F)
write.table(gene_mRNA_nodup, "snapfun_gene_mRNA_nodup.tsv", quote=F, col.names=T, row.names=T)
