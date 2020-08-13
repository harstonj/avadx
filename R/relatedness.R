#!/usr/bin/env Rscript

library(optparse, quietly=T, warn.conflicts=F)
library(SNPRelate, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)

option_list = list(
  make_option(c("-i", "--vcf_file_path"), type="character", default=NULL, 
              help="the path to the input vcf file; *.vcf.gz or *.vcf", metavar="character"),
  make_option(c("-g", "--gds_file_path"), type="character", default=NULL, 
              help="the path to the gds file; *.gds", metavar="character"),
  make_option(c("-c", "--kinship_cutoff"), type="character", default=0.3, 
              help="a numeric cutoff indicating 'related' by the kinship value; default is 0.3", metavar="character"),
  make_option(c("-o", "--output_folder_path"), type="character", default=NULL, 
              help="the path to the output folder", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="the path to the output folder", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

f = opt$vcf_file_path
g = opt$gds_file_path
threads = opt$threads
snpgdsVCF2GDS(f, g, method="biallelic.only")

genofile <- snpgdsOpen(g)

# Prune LD:
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, num.thread=threads)
# Get all selected snp id
snpset.id <- unlist(snpset)

# IBD (identity-by-descent):
ibd <- snpgdsIBDMoM(genofile, snp.id=snpset.id, maf=.05, missing.rate=.05, num.thread=threads)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

p <- ggplot(ibd.coeff, aes(x=kinship)) + geom_histogram(fill="yellow") + stat_bin(binwidth=.02, fill="#3399CC") + ylim(c(0, nrow(ibd.coeff))) + stat_bin(binwidth=.02, geom="text", aes(label=..count..), vjust=-2) + ylab("Number of individual pairs") + theme_classic(base_size=15)

setwd(opt$output_folder_path)
ggsave("IBD_histogram.pdf", p, width=9, height=9)

write.table(ibd.coeff, "IBD.txt", quote=F, sep="\t", row.names=F)
snpgdsClose(genofile)

related_pairs <- ibd.coeff[ibd.coeff$kinship > opt$kinship_cutoff,]
write.table(related_pairs, "IBD_related.txt", quote=F, sep="\t", row.names=F)
