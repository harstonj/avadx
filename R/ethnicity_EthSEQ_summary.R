#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

report <- args[1]
sampleID <- args[2]
out.dir <- args[3]

df <- read.table(report, header=T, stringsAsFactors=F, sep="\t")
sample.ID <- read.table(sampleID, header=F, stringsAsFactors=F)
sample.ID <- sample.ID$V1

df$sample.id <- sample.ID
rm(sample.ID, sampleID)

print(table(df$pop, df$type))
num_eur_inside <- sum(df$pop=="EUR" & df$type=="INSIDE")
num_eur_closest <- sum(df$pop=="EUR" & df$type=="CLOSEST")
print(paste0("Number of individuals inside EUR annotation: ", num_eur_inside, " (", round(num_eur_inside/nrow(df)*100), "%)."))
print(paste0("Number of individuals closest EUR annotation: ", num_eur_closest, " (", round(num_eur_closest/nrow(df)*100), "%)."))

setwd(out.dir)
write.table(df$sample.id[df$pop=="EUR" & df$type=="INSIDE"], "sampleID_inside_EUR.txt", quote=F, col.names=F, row.names=F)
write.table(df$sample.id[df$pop=="EUR" & df$type=="CLOSEST"], "sampleID_closest_EUR.txt", quote=F, col.names=F, row.names=F)
