import sys
import gzip


vcf_input = sys.argv[1]
vcf_output = sys.argv[2]
cutoff_AB_low = float(sys.argv[3]) if len(sys.argv) >= 4 else 0.3   # genotype with (calculated) AB < cutoff_AB_low will be converted into './.'
cutoff_AB_high = float(sys.argv[4]) if len(sys.argv) >= 5 else 0.7  # genotype with (calculated) AB > cutoff_AB_high will be converted into './.'
cutoff_DP = int(sys.argv[5]) if len(sys.argv) >= 6 else 4           # genotype with DP < cutoff_DP will be converted into './.'
cutoff_GQ = int(sys.argv[6]) if len(sys.argv) >= 7 else 15          # genotype with GQ < cutoff_GQ will be converted into './.'
cutoff_MR = float(sys.argv[7]) if len(sys.argv) >= 8 else 0.20      # variant sites with Missing Rate
line_buffer = 1e4

with (lambda f: (gzip.open if f.endswith('.gz') else open)(f, 'rt'))(vcf_input) as vcf_handle, gzip.open(vcf_output, 'wt') as fw, open(vcf_output + '.log', 'w') as wlog:
    line_cnt = 0
    for line in vcf_handle:
        if line[0] != '#':
            line_cnt += 1
            row = line.strip().split('\t')
            fm = row[8]
            fm_list = fm.split(':')
            idv = row[9:]
            idv_list = [x.split(':') for x in idv]
            has_GT, has_AD, has_DP, has_GQ, has_PL = [_ in fm_list for _ in ['GT', 'AD', 'DP', 'GQ', 'PL']]

            GT_idx = fm_list.index('GT') if has_GT else -1
            gt = [x[GT_idx] for x in idv_list] if GT_idx >= 0 else []

            AD_idx = fm_list.index('AD') if has_AD else -1
            ad = [x[AD_idx] for x in idv_list] if AD_idx >= 0 else []

            DP_idx = fm_list.index('DP') if has_DP else -1
            dp = [-1 if x[DP_idx] == '.' else x[DP_idx] for x in idv_list] if DP_idx >= 0 else []  # -1 means no call

            GQ_idx = fm_list.index('GQ') if has_GQ else -1
            gq = [-1 if x[GQ_idx] == '.' else x[GQ_idx] for x in idv_list] if GQ_idx >= 0 else []   # -1 means no call

            PL_idx = fm_list.index('PL') if has_PL else -1
            pl = [x[PL_idx] for x in idv_list] if PL_idx >= 0 else []

            # ab is a list:
            #   0 means hom ref
            #   1 means hom alt
            #  -1 means no call
            #  -2 means multi-allelic
            #  -3 means error (sum of all ADs is 0)
            #  99 means hom ref
            ab = []

            # Calculate AB
            if len(ad) != 0:  # Calculate AB from AD
                ad = [x.split(',') for x in ad]
                idx = 0
                for a in ad:  # a: ['0', '37']
                    GT = gt[idx]
                    if '.' in a or './.' in a:
                        ab.append(-1)  # no call
                    else:
                        a = [int(x) for x in a]
                        if len(a) > 2:  # multi-allelic loci
                            ab.append(-2)  # multi-allelic
                        elif len(a) == 2 and sum(a) != 0:
                            if GT == '0/1':
                                ab_temp = float(a[1]) / float(sum(a))
                                ab.append(ab_temp)
                            else:
                                ab.append(99)  # homozygous calls
                        else:
                            ab.append(-3)  # sum(AD) == 0 or others
                    idx += 1
                # Count bad genotypes
                n_bad_GTs = 0
                for i in range(0, len(gt)):
                    check_AB = float(ab[i])
                    check_DP = float(dp[i])
                    check_GQ = float(gq[i])
                    if check_AB < cutoff_AB_low or (check_AB < 1 and check_AB > cutoff_AB_high) or check_DP < cutoff_DP or check_GQ < cutoff_GQ:
                        n_bad_GTs += 1
                mr = n_bad_GTs / float(len(gt))
                if mr < cutoff_MR:
                    fw.write(line)
                    if line_cnt % 1e4 == 0:
                        fw.flush()
                else:
                    wlog.write('\t'.join(row[0:5]) + '\tHigh_Missing_Rate\t' + str(mr) + '\n')
            elif len(ad) == 0:  # No AB and no AD available.
                wlog.write('\t'.join(row[0:5]) + '\tNo_AB_or_AD_field\t' + 'NA' + '\n')
        else:
            fw.write(line)
