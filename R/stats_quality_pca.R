#!/usr/bin/env Rscript

library(optparse, quietly=T, warn.conflicts=F)
library(data.table, quietly=T, warn.conflicts=F)
library(ggfortify, quietly=T, warn.conflicts=F)

option_list = list(
  make_option(c("-f", "--stats"), type="character", default=NULL, 
              help="the path to the stats file output from bcftools", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="path for output pdf", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

f <- opt$stats

system(paste0("grep PSC ", f," > ", f, ".PSC.txt"))
df <- fread(paste0(f, ".PSC.txt"), skip=1, data.table=F)

rownames(df) <- df$`[3]sample`
df_ = df[, c(4:8,10,11,14)]
df_ = df_[ , which(apply(df_, 2, var) != 0)]
pca <- prcomp(df_, scale.=T)

p <- autoplot(pca, label=T, label.size=4, loadings=T, loadings.label=T, label.vjust=1.5, loadings.label.vjust=1.5)

ggsave(filename=opt$output,
	   plot=p,
	   width=9, height=9)
