import sys
import gzip
import re
# Yanran Wang
# Use:
#    python filterVCF_by_ABAD.py input.vcf.gz output.vcf.gz

cutoff_AB_low = 0.3   # genotype with (calculated) AB < cutoff_AB_low will be converted into "./."
cutoff_AB_high = 0.7  # genotype with (calculated) AB > cutoff_AB_high will be converted into "./."
cutoff_DP = 4          # genotype with DP < cutoff_DP will be converted into "./."
cutoff_GQ = 15         # genotype with GQ < cutoff_GQ will be converted into "./."
cutoff_MR = 0.20        # variant sites with Missing Rate

xopen = lambda f: (gzip.open if f.endswith(".gz") else open)(f)
vcf = xopen(sys.argv[1])

fw_p = sys.argv[2]
fw = gzip.open(fw_p, "wb")

flog_p = fw_p + ".log"
wlog = open(flog_p, "w")

c = 0
for line in vcf:
	line = line.decode("utf-8") if type(line) == bytes else line
	if line[0] != "#":
		l = (line.strip()).split("\t")
		position = "\t".join(l[0:5])
		#print position

		fm = l[8]
		fm_list = fm.split(":")

		idv = l[9:]
		idv_list = [x.split(":") for x in idv]

		GT_idx = fm_list.index("GT")
		gt = [x[GT_idx] for x in idv_list]
		try:
			AD_idx = fm_list.index("AD")
			ad = [x[AD_idx] for x in idv_list]
			#print 'ad', ad
		except ValueError:
			AD_idx = -1
			ad = []
		try:
			DP_idx = fm_list.index("DP")
			dp = [x[DP_idx] for x in idv_list]
			dp = [x if x!="." else -1 for x in dp] # -1 means no call
			#print 'dp', dp
		except ValueError:
			DP_idx = -1
			dp = []
		try:
			GQ_idx = fm_list.index("GQ")
			gq = [x[GQ_idx] for x in idv_list]
			gq = [x if x!="." else -1 for x in gq] # -1 means no call
			#print "gq", gq
		except ValueError:
			GQ_idx = -1
			gq = []
		try:
			PL_idx = fm_list.index("PL")
			pl = [x[PL_idx] for x in idv_list]
			#print "pl", pl
		except ValueError:
			PL_idx = -1
			pl = []
		#try:
		#	AB_idx = fm_list.index("AB")
		#	ab = [x[AB_idx] for x in idv_list]
		#	#print "ab", ab
		#except ValueError:
		#	AB_idx = -1
		#	ab = []

		ab = []
		# Calculate AB:
		if len(ad) != 0: # Calculate AB from AD
			#print position
			ad = [x.split(",") for x in ad]
			#print ad
			idx = 0
			for a in ad: # a: ['0', '37']
				GT = gt[idx]
				if "." in a or "./." in a:
					ab.append(-1) # no call
				else:
					a = [int(x) for x in a]
					if len(a) > 2: # multi-allelic loci
						ab.append(-2) # multi-allelic
					elif len(a)==2 and sum(a) != 0:
						if GT == "0/1":
							ab_temp = float(a[1]) / float(sum(a))
							ab.append(ab_temp)
						else:
							ab.append(99) # homozygous calls
					else:
						ab.append(-3) #sum(AD) == 0 or others
				idx += 1
				#print ab
				# ab is a list:
				#       0 means hom ref
				#       1 means hom alt
				#      -1 means no call
				#      -2 means multi-allelic
				#      -3 means error (sum of all ADs is 0)
				#      99 means hom ref

			# Count bad genotypes:
			n_bad_GTs = 0
			for i in range(0, len(gt)):
				check_AB = float(ab[i])
				check_DP = float(dp[i])
				check_GQ = float(gq[i])

				#print "AB ", check_AB
				#print "DP ", check_DP
				#print "GQ ", check_GQ

				if check_AB<cutoff_AB_low or (check_AB<1 and check_AB>cutoff_AB_high) or check_DP<cutoff_DP or check_GQ<cutoff_GQ:
					n_bad_GTs += 1

			mr = n_bad_GTs / float(len(gt))
			if mr < cutoff_MR:
				line = line if type(line) == bytes else line.encode("utf-8")
				fw.write(line)
			else:
				wlog.write(position + "\tHigh_Missing_Rate\t" + str(mr) + "\n")

		elif len(ad) == 0: # No AB and no AD available.
			wlog.write(position + "\tNo_AB_or_AD_field\t" + "NA" + "\n")
			#print "Cannot calculate AB!"
		#print ab

		c += 1

		#print c
		#if c == 100: break
	else:
		line = line if type(line) == bytes else line.encode("utf-8")
		fw.write(line)

vcf.close()
fw.close()
wlog.close()
