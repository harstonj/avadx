#/bin/bash

exome_file=$1
genome_file=$2
data_folder=$3

gunzip ${exome_file}.idx.gz
rm -rf ${data_folder}/gnomad_exome_allAFabove0.tmp_splits && mkdir ${data_folder}/gnomad_exome_allAFabove0.tmp_splits
cd ${data_folder}/gnomad_exome_allAFabove0.tmp_splits
zcat ${exome_file}.gz | split -d --lines=50000000

gunzip ${genome_file}.idx.gz
rm -rf ${data_folder}/gnomad_genome_allAFabove0.tmp_splits && mkdir ${data_folder}/gnomad_genome_allAFabove0.tmp_splits
cd ${data_folder}/gnomad_genome_allAFabove0.tmp_splits
zcat ${genome_file}.gz | split -d --lines=50000000
