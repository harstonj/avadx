#!/usr/bin/env Rscript

## This script is for merging all .gs files into one dataframe.
## Input of this script should be the folder name, which contains .gs files for all samples
## This script automatically check if the input file contains "normalized" or "un-normalized" gene scores from each sample file and outputs corresponding tables
## The output table contains these columns: Gene name, Transcript ID, Sample 1, Sample 2, ..., Sample N.
## There is NO label column

library(optparse, quietly=T, warn.conflicts=F)
library(data.table, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(purrr, quietly=T, warn.conflicts=F)

option_list = list(
  make_option(c("-f", "--input_file_folder"), type="character", default=NULL, 
              help="the path to the folder which contains all .gs files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="path to the output folder", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

# Read-in all .gs files into R:
setwd(opt$input_file_folder)

gs_list <- list.files(pattern=".gs")
gs <- list()
for(f in gs_list){
  gs[[gsub(".gs", "", f)]] <- fread(f, data.table=F)
  rm(f)
}


# un-normalized:
if("gene_score" %in% colnames(gs[[1]])){
  gs_unnormalized <- lapply(gs, function(x) {
    select(x,
           Gene, Transcript, gene_score)
  })
  gs_df_unnormalized <- gs_unnormalized %>% reduce(full_join, by=c("Gene", "Transcript"))
  colnames(gs_df_unnormalized)[3:ncol(gs_df_unnormalized)] <- names(gs_unnormalized)
  rm(gs_unnormalized)
  
  # Write to file:
  setwd(opt$out)
  fwrite(gs_df_unnormalized, "GeneScoreTable_raw.csv")
}


# normalized:
if("gene_score_normalized" %in% colnames(gs[[1]])){
  gs_normalized <- lapply(gs, function(x) {
    select(x,
           Gene, Transcript, gene_score_normalized)
  })
  gs_df_normalized <- gs_normalized %>% reduce(full_join, by=c("Gene", "Transcript"))
  colnames(gs_df_normalized)[3:ncol(gs_df_normalized)] <- names(gs_normalized)
  rm(gs_normalized)
  
  # Write to file
  setwd(opt$out)
  fwrite(gs_df_normalized, "GeneScoreTable_normalized.csv")
}

rm(gs)

