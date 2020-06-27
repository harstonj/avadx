#/bin/bash

data_folder=$1
exome_file=$2
genome_file=$3

rm -f ${data_folder}/gnomad_exome_allAFabove0.tmp
rm -f ${data_folder}/gnomad_genome_allAFabove0.tmp
rm -rf ${data_folder}/gnomad_genome_allAFabove0.tmp_splits

bgzip -f ${exome_file}
tabix -f ${exome_file}.gz

bgzip -f ${genome_file}
tabix -f ${genome_file}.gz
