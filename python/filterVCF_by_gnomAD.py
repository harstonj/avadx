import sys
import gzip
import subprocess
import csv
from pathlib import Path


vcf_input = sys.argv[1]
vcf_output = sys.argv[2]
gnomad_exome = sys.argv[3]
gnomad_genome = sys.argv[4]

base_out = Path(vcf_output).parent
vcf_regions = base_out / 'vcf_regions.tsv'
exome_filtered_path = base_out / 'gnomad_exome_filtered.tsv'
genome_filtered_path = base_out / 'gnomad_genome_filtered.tsv'

exome_filtered_dict = {}
genome_filtered_dict = {}

subprocess.call(f'zcat {vcf_input} | egrep -v "(^#.*|^$)" | cut -f1,2 > {vcf_regions}', shell=True)

subprocess.call(f'tabix {gnomad_exome} -R {vcf_regions} > {exome_filtered_path}', shell=True)
with exome_filtered_path.open() as fin:
    reader = csv.reader(fin, delimiter='\t')
    for row in reader:
        if not row[0] in exome_filtered_dict:
            exome_filtered_dict[row[0]] = set()
        exome_filtered_dict[row[0]].add(int(row[1]))

subprocess.call(f'tabix {gnomad_genome} -R {vcf_regions} > {genome_filtered_path}', shell=True)
with genome_filtered_path.open() as fin:
    reader = csv.reader(fin, delimiter='\t')
    for row in reader:
        if not row[0] in genome_filtered_dict:
            genome_filtered_dict[row[0]] = set()
        genome_filtered_dict[row[0]].add(int(row[1]))

vcf_handle = (lambda f: (gzip.open if f.endswith(".gz") else open)(f))(vcf_input)
vcf_filtered = gzip.open(vcf_output, "wb")

for line in vcf_handle:
    line_decoded = line.decode("utf-8")
    if line_decoded[0] != "#":
        chrom, pos = line_decoded.strip().split('\t', 2)[0:2]
        pos = int(pos)
        if pos in exome_filtered_dict.get(chrom, set()) or pos in genome_filtered_dict.get(chrom, set()):
            vcf_filtered.write(line)
        else:
            pass  # removed entries are not saved due to perfrmance considerations
    else:
        vcf_filtered.write(line)

vcf_handle.close()
vcf_filtered.close()
