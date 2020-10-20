#!/usr/bin/env Rscript
library(dplyr)
library(seqinr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

exonic.file <- args[1]
dbfolder <- args[2]
outputfolder <- args[3]

logs = c()

setwd(outputfolder)
all.file = gsub(".exonic_variant_function", ".variant_function", exonic.file)
all.info = data.table(table(fread(all.file, header=F, data.table=F)$V1))
colnames(all.info) = c('type', 'count')
all.info = rbindlist(list(all.info, data.table(type='TOTAL', count=sum(all.info$count))))
fwrite(all.info, 'annovar_stats.csv')

setwd(dbfolder)
df_len <- read.csv("Transcript-ProtLength.csv", stringsAsFactors=F, header=T)
#mutOut <- read.table("Mutations.mutOut", stringsAsFactors=F, header=F, sep="\t")

var_filter = c("synonymous SNV", "nonsynonymous SNV", "stopgain", "stoploss")
exonic <- fread(exonic.file, header=F, data.table=F)
logs = c(logs, paste('Exonic variants read:', nrow(exonic)))
exonic <- exonic[exonic$V3!="UNKNOWN",]
exonic <- exonic[exonic$V2 %in% var_filter, ]
vars_retained = table(exonic['V2'])
logs = c(logs, paste('Annotated variants retained:', nrow(exonic)), paste(names(vars_retained), collapse = ","), toString(vars_retained))

variants <- strsplit(exonic$V3, ",")
variants <- lapply(variants, function(x) strsplit(x, ":"))

variant_lengths <- lapply(variants, function(x) unlist(lapply(x, length)))
if(unique(unlist(variant_lengths))!=5){
  stop("ERROR: Annovar .exonic_variant_function file format not recognized.")
}

variants_mRNA <- lapply(variants, function(x) lapply(x, function(y) y[2]))
variants_aaMut <- lapply(variants, function(x) lapply(x, function(y) y[5]))
variants_gene <- lapply(variants, function(x) lapply(x, function(y) y[1]))

gene_mRNA_list <- list()
for(i in 1:length(variants)){
  g <- unlist(variants_gene[[i]])
  line <- exonic$V1[i]
  type <- exonic$V2[i]
  trans <- unlist(variants_mRNA[[i]])
  aaMut <- unlist(variants_aaMut[[i]])
  dfi <- data.frame("Gene"=g, "Transcript"=trans, "aaMut"=gsub("p.", "", aaMut), "line"=line, "type"=type)
  gene_mRNA_list[[i]] <- dfi
  rm(g, trans, dfi)
}

gene_mRNA <- rbindlist(gene_mRNA_list, use.names=T)
gene_mRNA$prot_length <- df_len$Prot_length[match(gene_mRNA$Transcript, df_len$Transcript)]
logs = c(logs, paste0("Variants across all isoforms: ", nrow(gene_mRNA), " (in ", length(unique(gene_mRNA$Transcript)), " transcripts)"))

setwd(outputfolder)
  write.table(unique(as.character(gene_mRNA$Transcript[is.na(gene_mRNA$prot_length)])), 
              file="missing_mRNA_accession_all_transcripts.txt", quote=F, col.names=F, row.names=F)
if(any(is.na(gene_mRNA$prot_length))){
  missing_count = length(which(is.na(gene_mRNA$prot_length)))
  missing_vars = gene_mRNA[which(is.na(gene_mRNA$prot_length)), ]
  missing_vars_cnt = length(unique(missing_vars$Transcript))
  logs = c(logs,
    paste0("Information missing for ", missing_vars_cnt, " transcripts (", round(missing_count/nrow(gene_mRNA)*100, 2), "% of all variants)")
  )
}

# check if Mutations.mutOut already contains the variants:
#gene_mRNA$Exist <- ifelse(paste0(gene_mRNA$Transcript, "-", gene_mRNA$aaMut) %in% paste0(mutOut$V1, "-", mutOut$V2), T, F)
#print(paste0(sum(gene_mRNA$Exist), " mutations (", round(sum(gene_mRNA$Exist)/nrow(gene_mRNA)*100),"%) already have SNAP scores in Mutations.mutOut. Another ", sum(!gene_mRNA$Exist), " new mutations need to run with SNAP."))
#print("Generating SNAP input files...")

# For the same line, keep only the mutation in the longest protein:
# -> if gene has multiple transcripts of identical lengths, keep the one with "smaller" NM accesion number as ordered by dplyr->arrange
#gene_mRNA <- gene_mRNA[gene_mRNA$Exist==F, ]
gene_mRNA_nodup <- gene_mRNA %>%
  group_by(line) %>%
  filter(prot_length == max(prot_length)) %>%
  arrange(Gene, Transcript) %>%
  filter(row_number()==1)
logs = c(logs, paste("\nFinal variants retained (longest isoform; primary annovar isoform if same length):", nrow(gene_mRNA_nodup)))

setwd(outputfolder)
write.table(unique(as.character(gene_mRNA_nodup$Transcript[is.na(gene_mRNA_nodup$prot_length)])), 
            file="missing_mRNA_accession_mapped_transcripts.txt", quote=F, col.names=F, row.names=F)
if(any(is.na(gene_mRNA_nodup$prot_length))){
  missing_count = length(which(is.na(gene_mRNA_nodup$prot_length)))
  missing_vars = gene_mRNA_nodup[which(is.na(gene_mRNA_nodup$prot_length)), ]
  missing_vars_cnt = length(unique(missing_vars$Transcript))
  logs = c(logs,
    paste0(
      "Information missing for ", missing_vars_cnt, " mapped transcripts (", round(missing_count/nrow(gene_mRNA_nodup)*100, 2), "% of all variants)"
    )
  )
}

# Read in fasta files:
setwd(dbfolder)
fas <- read.fasta("prot_seqs.fa", as.string=T, seqtype="AA", set.attributes=F)
fas_names <- names(fas)
# fas_mRNA <- strsplit(fas_names, "\\.|\\|")
# fas_mRNA <- unlist(lapply(fas_mRNA, function(x) x[2]))
fas_mRNA <- strsplit(fas_names, "\\.")
fas_mRNA <- unlist(lapply(fas_mRNA, function(x) x[1]))

# 1. check if the prot_seq.txt contains the protein sequence:
missing_prot_seq <- setdiff(unique(as.character(gene_mRNA$Transcript)), fas_mRNA)

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
  write.table(missing_prot_seq_nodup, "TranscriptAccess_missing_prot_seq_mapped.txt", quote=F, col.names=F, row.names=F)

if(length(missing_prot_seq_nodup)>0){
  print(paste0(length(missing_prot_seq_nodup), " QUERY transcripts do not have protein sequences in the db folder"))
}

# 2. Generate varidb input
setwd(outputfolder)

gene_mRNA_nodup_nonsyn = gene_mRNA_nodup[which(gene_mRNA_nodup$type == 'nonsynonymous SNV'), ]
matched_nonsyn = match(unique(gene_mRNA_nodup_nonsyn$Transcript), fas_mRNA)
write.fasta(fas[matched_nonsyn], fas_mRNA[matched_nonsyn], file.out='varidb_query_nonsyn.fa')
write.table(gene_mRNA_nodup_nonsyn[c('Transcript', 'aaMut')], "varidb_query_nonsyn.ids", quote=F, col.names=F, row.names=F)

gene_mRNA_nodup_syn = gene_mRNA_nodup[which(gene_mRNA_nodup$type == 'synonymous SNV'), ]
matched_syn = match(unique(gene_mRNA_nodup_syn$Transcript), fas_mRNA)
write.fasta(fas[matched_syn], fas_mRNA[matched_syn], file.out='varidb_query_syn.fa')
write.table(gene_mRNA_nodup_syn[c('Transcript', 'aaMut')], "varidb_query_syn.ids", quote=F, col.names=F, row.names=F)

gene_mRNA_nodup_stopgainloss = gene_mRNA_nodup[which(gene_mRNA_nodup$type %in% c('stopgain', 'stoploss')), ]
matched_stopgainloss = match(unique(gene_mRNA_nodup_stopgainloss$Transcript), fas_mRNA)
write.fasta(fas[matched_stopgainloss], fas_mRNA[matched_stopgainloss], file.out='varidb_query_stopgainloss.fa')
write.table(gene_mRNA_nodup_stopgainloss[c('Transcript', 'aaMut')], "varidb_query_stopgainloss.ids", quote=F, col.names=F, row.names=F)

write.table(gene_mRNA_nodup, "varidb_gene_mRNA_all.tsv", quote=F, col.names=T, row.names=T)

logs = c(logs, "\n-- Variant types --")
logs = c(logs, paste("synonymous SNVs    :", nrow(gene_mRNA_nodup_syn)))
logs = c(logs, paste("non-synonymous SNVs:", nrow(gene_mRNA_nodup_nonsyn)))
logs = c(logs, paste("stopgain / stoploss:", nrow(gene_mRNA_nodup_stopgainloss)))

if(length(missing_prot_seq)>0){
  setwd(outputfolder)
  write.table(missing_prot_seq, "TranscriptAccess_missing_prot_seq_all.txt", quote=F, col.names=F, row.names=F)
  logs = c(logs,
    paste0(
      "\n\n\nNOTE: For ", length(missing_prot_seq), " transcripts no protein sequences were available. ",
      "To include respective variants, add missing sequences to the avadx/prot_seqs.fa file. ",
      "mRNA transcript accession numbers can be found in the working directory: varidb/TranscriptAccess_missing_prot_seq_[all|mapped].txt. ",
      "Please use NCBI batch entrez query to obtain protein sequences (https://www.ncbi.nlm.nih.gov/sites/batchentrez)."
    )
  )
}

fileConn<-file("SNAP_extract_vars.log")
writeLines(logs, fileConn)
close(fileConn)
