library(data.table)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

data.folder = args[1]
hgref = args[2]

exome.file = paste0(data.folder, '/gnomad_exome_allAFabove0.tmp')
genome.file = paste0(data.folder, '/gnomad_genome_allAFabove0.tmp')

exome.file.out = paste0(data.folder, "/", hgref, "_gnomad_exome_allAFabove0.txt")
gnomad.exome_dt = data.table::fread(exome.file, header=T, data.table=T)
gnomad.exome_dt_filtered = gnomad.exome_dt[gnomad.exome_dt$gnomAD_exome_ALL > 0, c('#Chr', 'Start')]
rm(gnomad.exome_dt)
data.table::fwrite(gnomad.exome_dt_filtered, file=exome.file.out,  col.names=FALSE, sep="\t")
rm(gnomad.exome_dt_filtered)

genome.file.out = paste0(data.folder, "/", hgref, "_gnomad_genome_allAFabove0.txt")
gnomad.genome.splits = list.files(path = paste0(data.folder, "/gnomad_genome_allAFabove0.tmp_splits"), full.names=T)

has_header=T
for (genome.file.split in gnomad.genome.splits) {
    print(genome.file.split)
    gnomad.genome_dt = data.table::fread(genome.file.split, header=has_header, data.table=T, showProgress=F)
    gnomad.genome_dt_filtered = gnomad.genome_dt[gnomad.genome_dt[[3]] > 0, 1:2]
    rm(gnomad.genome_dt)
    data.table::fwrite(gnomad.genome_dt_filtered, file=genome.file.out, col.names=FALSE, append=!has_header, sep="\t")
    rm(gnomad.genome_dt_filtered)
    has_header=F
}
