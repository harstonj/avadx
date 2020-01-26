#!/usr/bin/env Rscript
library(EthSEQ)

args <- commandArgs(trailingOnly=TRUE)

# Sys.setenv('R_MAX_VSIZE'=32000000000)

ethseq.Analysis(
  target.vcf = args[1],
  out.dir = args[2],
  model.available = "SS2.Major",
  verbose = TRUE,
  composite.model.call.rate = 1,
  space = "3D"
)
