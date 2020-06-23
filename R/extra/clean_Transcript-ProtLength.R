#!/usr/bin/env Rscript

# This script is used to clean the Transcript-ProtLength.csv file in the db folder
# The purpose of the clean is to keep only the transcript with the longest protein
# If one gene has multiple transcript with the same protein lengths, we arbitrarily choose to keep the transcript with "smaller" NM accesion number by dplyr arrange

library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

file_path <- args[1]
output_path <- args[2]

# Read-in Acc-protein_length:
df <- fread(file_path, data.table=F)
#   keep the transcript with longer protein for the same gene:
df_flted <- df %>% 
  group_by(Gene) %>% top_n(1, Prot_length) %>%
  arrange(Gene, Transcript) %>%
  group_by(Gene) %>% filter(row_number()==1)

if(any(duplicated(df_flted$Gene))==F){
  write.csv(df_flted, output_path, quote=F, row.names=F)
}

