from Bio import SeqIO
import csv

seq_dict = {}
for seq_record in SeqIO.parse("/Users/WangYanran/Documents/BrombergLab/SNAP_database/NCBI/GCF_000001405.25_GRCh37.p13_rna.gbff", "genbank"):
	#print(repr(seq_record.seq))
	#print(len(seq_record))
	
	for feature in seq_record.features:
		if 'protein_id' in feature.qualifiers:
			#print(seq_record.id)
			#print(feature.qualifiers['protein_id'][0])
			#print("\n")
			seq_dict[seq_record.id] = feature.qualifiers['protein_id'][0]

with open("/Users/WangYanran/Documents/BrombergLab/SNAP_database/parse_result.csv", "w") as f:
	for key in seq_dict.keys():
		f.write("%s,%s\n"%(key, seq_dict[key]))
