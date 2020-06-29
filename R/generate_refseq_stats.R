library(data.table)
library(tidyr)
library(dplyr)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

stats.file = args[1]
seq.file = args[2]
data.folder = args[3]

stats_dt = data.table::fread(stats.file, header=T, data.table=T, na.strings="")
stats_dt = stats_dt[stats_dt$'# feature' == 'CDS' & stats_dt$class == 'with_protein'][, c('symbol', 'product_accession', 'related_accession', 'feature_interval_length', 'product_length', 'name')]
stats_dt_filtered = stats_dt[!is.na(stats_dt$related_accession), ]
stats_dt_sub = stats_dt_filtered[, c('symbol', 'related_accession', 'product_length')]
setnames(stats_dt_sub, 1:3, c('Gene', 'Transcript', 'Prot_length'))
stats_dt_sub[, Transcript := tstrsplit(Transcript, ".", fixed=TRUE, keep=c(1))]
fwrite(stats_dt_sub, file = paste0(data.folder, "/Transcript-ProtLength.csv"))
mapping_dt = stats_dt_filtered[, c('product_accession', 'related_accession')]
setnames(mapping_dt, 1:2, c('protein', 'transcript'))
fwrite(mapping_dt, file = paste0(data.folder, "/refseq_mapping.csv"))

seq.fasta = read.fasta(seq.file, as.string=T, seqtype="AA", set.attributes=F)
seq.names <- names(seq.fasta)
matched = match(seq.names, mapping_dt$protein)
write.fasta(seq.fasta[!is.na(matched)], mapping_dt$transcript[matched], file.out=paste0(data.folder, '/prot_seqs.fa'))

# Clean the Transcript-ProtLength.csv file -  keep only the transcript with the longest protein
# If one gene has multiple transcript with the same protein lengths,
# we arbitrarily choose to keep the transcript with "smaller" NM accesion number by dplyr arrange

#   keep the transcript with longer protein for the same gene:
df_flted <- stats_dt_sub %>% 
  group_by(Gene) %>% top_n(1, Prot_length) %>%
  arrange(Gene, Transcript) %>%
  group_by(Gene) %>% filter(row_number()==1)

if(any(duplicated(df_flted$Gene))==F){
  fwrite(df_flted, file = paste0(data.folder, "/Transcript-ProtLength.csv"))
}
