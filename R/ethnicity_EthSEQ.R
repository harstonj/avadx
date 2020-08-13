#!/usr/bin/env Rscript
library(EthSEQ)

args <- commandArgs(trailingOnly=TRUE)

# Sys.setenv('R_MAX_VSIZE'=32000000000)

ethseq.Analysis(
  target.vcf = args[1],
  model.folder = args[2],
  out.dir = args[3],
  model.available = "Exonic.All",
  verbose = TRUE,
  composite.model.call.rate = 1,
  space = "3D",
  cores = strtoi(args[4])
)
