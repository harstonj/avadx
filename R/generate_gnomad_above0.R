library(dplyr)
library(data.table)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

data.folder = args[1]
hgref = args[2]
chr_order = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

for (type in c('exome', 'genome')) {
    file.out = paste0(data.folder, "/", hgref, "_gnomad_", type, "_allAFabove0.txt")
    gnomad.splits = list.files(path = paste0(data.folder, "/gnomad_", type, "_allAFabove0.tmp_splits"), full.names=T)
    first_split=T
    add_from_next = 0
    for (i in seq(1, length(gnomad.splits))) {
        # process current split
        file.split = gnomad.splits[i]
        print(file.split)
        gnomad_dt = data.table::fread(file.split, header=first_split, data.table=T, showProgress=F, stringsAsFactors=T, select=list(factor=1, integer=2, numeric=6), na.strings=".", skip=add_from_next)
        gnomad_dt_filtered = gnomad_dt[gnomad_dt[[3]] > 0, 1:2]
        rm(gnomad_dt)
        levels_order = intersect(chr_order, levels(gnomad_dt_filtered[[1]]))
        chr_col = colnames(gnomad_dt_filtered)[1]
        gnomad_dt_filtered[, (chr_col) := factor(get(chr_col), levels = levels_order)]
        gnomad_dt_filtered_ordered = setDT(gnomad_dt_filtered %>% group_by_at(1) %>% arrange_at(2, .by_group=T))
        rm(gnomad_dt_filtered)
               
        # check next split
        last_chr = paste(gnomad_dt_filtered_ordered[[nrow(gnomad_dt_filtered_ordered), 1]])
        last_pos = as.integer(gnomad_dt_filtered_ordered[[nrow(gnomad_dt_filtered_ordered), 2]])
        if (!first_split & i+1 <= length(gnomad.splits)) {
            following_lines = scan(gnomad.splits[i+1], nlines = 2, what = character(), sep = "\n")
            add_from_next = 0
            for (line_i in following_lines) {
                line_tsv = strsplit(line_i, "\t")[[1]]
                line_chr = paste(line_tsv[1])
                line_pos = as.integer(line_tsv[2])
                if (line_chr != last_chr) {
                    next
                } else if (line_pos < last_pos) {
                    add_from_next = add_from_next + 1
                }
            }
            if (add_from_next > 0) {
                last_current = gnomad_dt_filtered_ordered[nrow(gnomad_dt_filtered_ordered), ]
                gnomad_dt_filtered_ordered = gnomad_dt_filtered_ordered[-nrow(gnomad_dt_filtered_ordered), ]
                data.table::fwrite(gnomad_dt_filtered_ordered, file=file.out, col.names=first_split, append=!first_split, sep="\t")
                for (add_i in seq(1,add_from_next)) {
                    line_tsv = strsplit(following_lines[add_i], "\t")[[1]]
                    if (as.numeric(line_tsv[6]) > 0) {
                        add_dt = data.table(line_tsv[1], line_tsv[2])
                        names(add_dt) = names(gnomad_dt_filtered_ordered)
                        data.table::fwrite(add_dt, file=file.out, col.names=F, append=T, sep="\t")
                    }
                }
            } else {
                data.table::fwrite(gnomad_dt_filtered_ordered, file=file.out, col.names=first_split, append=!first_split, sep="\t")
            }
        } else {
            data.table::fwrite(gnomad_dt_filtered_ordered, file=file.out, col.names=first_split, append=!first_split, sep="\t")
        }
        rm(gnomad_dt_filtered_ordered)
        first_split=F
    }
}
