#/bin/bash

exome_file=$1
genome_file=$2
data_folder=$3

rm -rf ${data_folder}/gnomad_exome_allAFabove0.tmp_splits && mkdir ${data_folder}/gnomad_exome_allAFabove0.tmp_splits
split -d --lines=50000000 ${exome_file} ${data_folder}/gnomad_exome_allAFabove0.tmp_splits/

rm -rf ${data_folder}/gnomad_genome_allAFabove0.tmp_splits && mkdir ${data_folder}/gnomad_genome_allAFabove0.tmp_splits
split -d --lines=50000000 ${genome_file} ${data_folder}/gnomad_genome_allAFabove0.tmp_splits/
