#!/usr/bin/env Rscript

## This script is for merging all .gs files into one dataframe.
## Input of this script should be the folder name, which contains .gs files for all samples
## This script automatically check if the input file contains "normed" or "un-normed" gene scores from each sample file and outputs corresponding tables
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

# Un-normed:
if("gene_score" %in% colnames(gs[[1]])){
  gs_unnormed <- lapply(gs, function(x) {
    select(x,
           Gene, Transcript, gene_score)
  })
  gs_df_unnormed <- gs_unnormed %>% reduce(full_join, by=c("Gene", "Transcript"))
  colnames(gs_df_unnormed)[3:ncol(gs_df_unnormed)] <- names(gs_unnormed)
  rm(gs_unnormed)
  
  # Write to file:
  setwd(opt$out)
  fwrite(gs_df_unnormed, "GeneScoreTable_unnormed.txt")
}


# Normed:
if("gene_score_normed" %in% colnames(gs[[1]])){
  gs_normed <- lapply(gs, function(x) {
    select(x,
           Gene, Transcript, gene_score_normed)
  })
  gs_df_normed <- gs_normed %>% reduce(full_join, by=c("Gene", "Transcript"))
  colnames(gs_df_normed)[3:ncol(gs_df_normed)] <- names(gs_normed)
  rm(gs_normed)
  
  # Write to file
  setwd(opt$out)
  fwrite(gs_df_normed, "GeneScoreTable_normed.txt")
}

rm(gs)

