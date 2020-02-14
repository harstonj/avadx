#!/usr/bin/env Rscript

## This script is for gene_score calculation.
## Input and output of this script is only for one .exonic_variant_function file, i.e. one individual/sample
## For calculation of multiple samples, user can either:
##  (1) Use a for loop in the shell script which calls this Rscript.
##  (2) Use cluster (amarel) and run this paralelly for each input files using arrays.

library(optparse, quietly=T, warn.conflicts=F)
library(data.table)
library(tidyr)
library(dplyr)

option_list = list(
  make_option(c("-f", "--exonic_file"), type="character", default=NULL, 
              help=".exonic_variant_function by ANNOVAR of one individual", metavar="character"),
  make_option(c("-s", "--snapscore_file"), type="character", default=NULL, 
              help="snapscore file (Mutations.mutOut) where columns are: gene (NM_ number), aa_mutation, snapscore", metavar="character"),
  make_option(c("-l", "--protlength_file"), type="character", default=NULL, 
              help="protein length file where columns are: gene name, NM_ number, corresponding protein length from RefSeq", metavar="character"),
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="method to use for calculating gene score: 'sum' - as the prototype; 'production' - as in the VTE paper", metavar="character"),
  make_option(c("-n", "--norm"), type="character", default=NULL, 
              help="whether to normalize the gene score by protein length: 'y' - yes; 'n' - no; 'both' - return both results of yes and no. Note that the normalization is done by: gene score * 100 / protein length, so that the normalized score ranges approximately from 0 to 1.", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="path to the output folder. The script outputs two files: .txt file for a summary of counts; .gs file for the gene score table", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

# Extract input file name:
sample_name <- gsub(".avinput.exonic_variant_function", "", basename(opt$exonic_file))

# Read-in SNAP scores:
snap_score <- fread(opt$snapscore_file, data.table=F)
colnames(snap_score) <- c("Transcript", "Mut", "SNAP_Score")

# Read-in Acc-protein_length:
protL <- fread(opt$protlength_file, data.table=F)
protL_new <- protL %>% select(Transcript, Prot_length)

# Pre-precessing of the input file:
exonic <- fread(opt$exonic_file, data.table=F)
exonic_flted <- exonic %>% 
  filter(V3 != "UNKNOWN") %>%
  select(V1, V2, V3, V9) %>%
  separate_rows(V3, sep=",") %>%
  filter(nchar(V3)!=0) %>%
  separate(V3, c("Gene", "Transcript", "exon", "nuChange", "aaChange"), sep=":") %>%
  left_join(protL_new, by="Transcript") %>%
  filter(!is.na(Prot_length)) %>%
  mutate(Mut=gsub("p.", "", aaChange)) %>%
  left_join(snap_score, by=c("Transcript", "Mut")) %>%
  mutate(v_score = ifelse(is.na(SNAP_Score) & V2=="nonsynonymous SNV", NA, 
                          ifelse(is.na(SNAP_Score) & (V2=="stoploss" | V2=="stopgain"), 100, 
                                 ifelse(is.na(SNAP_Score) & V2=="synonymous SNV", -99, SNAP_Score)))) %>%
  filter(!is.na(v_score))

# Generate a summary for this sample:
exonic_flted_summary <- exonic_flted[,c("V2", "V9", "v_score")]
exonic_flted_summary <- exonic_flted_summary %>%
  mutate(neutrality = ifelse(v_score > 0 & V2=="nonsynonymous SNV", "NonNeutral", ifelse(v_score <= 0 & V2=="nonsynonymous SNV", "Neutral", NA)))

summary_exonic <- as.data.frame(table(exonic_flted_summary$V2, exonic_flted_summary$V9, exonic_flted_summary$neutrality, exclude=NULL))
summary_exonic <- summary_exonic %>% filter(Freq != 0)

setwd(opt$out)
write.table(summary_exonic, paste0("summary_", sample_name, ".txt"), sep="\t", quote=F, col.names=F, row.names=F)


# Score calculation:
if(opt$m=="sum"){
  
  # Generate variant scores:
  df <- exonic_flted %>%
    mutate(heti = ifelse(V9=="het", 0.25, 1)) %>%
    mutate(variant_score = ifelse(V2 == "synonymous SNV", 0.05 * heti,
                                  ifelse(V2 == "nonsynonymous SNV" & v_score <= 0, 0.055 * heti,
                                         ifelse(v_score > 0, (0.06+(v_score/100*.94)) * heti, 
                                                NA)))) 
  # Summing variant scores within each gene:
  df_gs <- df %>%
    group_by(Gene, Transcript) %>% summarise(gene_score = sum(variant_score))
  #    the larger the gene_score, the greater the function change
  
  # Normalizing variant score by protein lengths:
  df_gs$Prot_length <- protL$Prot_length[match(df_gs$Transcript, protL$Transcript)]
  df_gs <- df_gs %>%
    mutate(gene_score_normed = gene_score * 100/Prot_length)


  }else if(opt$m=="production"){
    # Generate variant scores:
    df <- exonic_flted %>%
      mutate(heti = ifelse(V9=="het", 0.25, 1)) %>%
      mutate(variant_score = ifelse(V2 == "synonymous SNV", 0.05 * heti,
                                    ifelse(V2 == "nonsynonymous SNV" & v_score <= 0, 0.055 * heti,
                                           ifelse(v_score > 0, (0.06+(v_score/100*.94)) * heti, 
                                                  NA)))) 
    # Summing variant scores within each gene:
    df_gs <- df %>%
      group_by(Gene, Transcript) %>% summarise(gene_score = prod((1 - variant_score))) 
    #    gene_score of 0 means protein function totally disrupted (e.g. has stopgain/loss)
    #    gene_score close to 1 means protein function completely kept (e.g. has only synonymous variant or neutral nonsynonymous variant)
    
    # Normalizing variant score by protein lengths:
    df_gs$Prot_length <- protL$Prot_length[match(df_gs$Transcript, protL$Transcript)]
    df_gs <- df_gs %>%
      mutate(gene_score_normed = gene_score * 100/Prot_length)  
}

setwd(opt$out)

if(opt$n=="n"){
  fwrite(df_gs[, c("Gene", "Transcript", "gene_score")], paste0(sample_name, ".gs"))
}else if(opt$n=="y"){
  fwrite(df_gs[, c("Gene", "Transcript", "gene_score_normed")], paste0(sample_name, ".gs"))
}else if(opt$n=="both"){
  fwrite(df_gs[, c("Gene", "Transcript", "gene_score", "gene_score_normed")], paste0(sample_name, ".gs"))
}


