#/bin/bash

exome_file=$1
genome_file=$2
data_folder=$3

cut -f 1,2,6 ${exome_file} > ${data_folder}/gnomad_exome_allAFabove0.tmp
cut -f 1,2,6 ${genome_file} > ${data_folder}/gnomad_genome_allAFabove0.tmp

rm -rf ${data_folder}/gnomad_genome_allAFabove0.tmp_splits && mkdir ${data_folder}/gnomad_genome_allAFabove0.tmp_splits
split -d --lines=50000000 ${data_folder}/gnomad_genome_allAFabove0.tmp ${data_folder}/gnomad_genome_allAFabove0.tmp_splits/
