#/bin/bash

exome_file=$1
genome_file=$2
data_folder=$3

rm -rf ${data_folder}/gnomad_exome_allAFabove0.tmp_splits
rm -rf ${data_folder}/gnomad_genome_allAFabove0.tmp_splits

bgzip -f ${exome_file}
tabix -f -b 2 ${exome_file}.gz

bgzip -f ${genome_file}
tabix -f -b 2 ${genome_file}.gz
