import sys
import gzip
#import re
import subprocess
# Yanran Wang
# Use:
#   python filterVCF_by_gnomAD.py input.vcf.gz output.vcf.gz

xopen = lambda f: (gzip.open if f.endswith(".gz") else open)(f)
vcf = xopen(sys.argv[1])

fw_p = sys.argv[2]
fw = gzip.open(fw_p, "wb")

gnomad_exome = sys.argv[3]
gnomad_genome = sys.argv[4]

flog_p = fw_p + ".removed.gz"
wlog = gzip.open(flog_p, "wb")

c = 0
for line in vcf:
	line = line.decode("utf-8") if type(line) == bytes else line
	if line[0] != "#":
		c += 1
		l = (line.strip()).split("\t")
		#position = "\t".join(l[0:5])
		chrom = l[0]
		pos = l[1]
		rs = l[2]
		ref = l[3]
		alt = l[4]
		#print position
		#fm = l[8]
		#fm_list = fm.split(":")
		#idv = l[9:]
		#idv_list = [x.split(":") for x in idv]
		search_res_exome = subprocess.check_output("tabix {} {}:{}-{}".format(gnomad_exome, chrom, pos, int(pos)+1), shell=True)
		search_res_genome = subprocess.check_output("tabix {} {}:{}-{}".format(gnomad_genome, chrom, pos, int(pos)+1), shell=True)
		if (pos in search_res_exome) or (pos in search_res_genome):
			#print "EXIST pos: {}-{}, ref: {}, alt: {}".format(chrom, pos, ref, alt)
			line = line if type(line) == bytes else line.encode("utf-8")
			fw.write(line)
		else:
			#print "NON-EXIST pos: {}-{}, ref: {}, alt: {}".format(chrom, pos, ref, alt)
			line = line if type(line) == bytes else line.encode("utf-8")
			wlog.write(line)

		if c%1000==0: print(c)
	else:
		line = line if type(line) == bytes else line.encode("utf-8")
		fw.write(line)

vcf.close()
fw.close()
wlog.close()
