#!usr/bin/python

# usage:
# python snap-scores.py snap-scores/RALB
# output appends to file "Mutations.mutOut"

import sys
from os.path import basename
fn = sys.argv[1]

gene = basename(fn)
gene = gene.replace(".out", "")

fw = open("Mutations.mutOut", "a")

with open(fn, "r") as f:
	for line in f:
		line = line.split("|")
		if len(line) > 1:
			#print line
			mut = line[0]
			score = line[len(line)-1]

			mut = mut[:mut.find("=>")]
			score = score[(score.find("=")+1):]

			mut = mut.strip()
			score = score.strip()
			#print gene, mut, score
			#print gene + "\t" + mut + "\t" + score + "\n"
			fw.write(gene + "\t" + mut + "\t" + score + "\n")

fw.close()

