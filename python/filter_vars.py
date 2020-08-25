import sys
import pandas as pd
from pathlib import Path
from avadx.helper import parse_fasta


SNAP_LIMIT = 6000

wd_path = Path(sys.argv[1])
non_syn_ids = wd_path / 'varidb_query_nonsyn.ids'
non_syn_fa = wd_path / 'varidb_query_nonsyn.fa'
non_syn_ids_filtered = wd_path / 'varidb_query_nonsyn_filtered.ids'
non_syn_fa_filtered = wd_path / 'varidb_query_nonsyn_filtered.fa'

report_name = 'varidb_prefilter_report.log'
excluded_name = 'varidb_prefilter_excluded.tsv'

non_syn_vars = pd.read_table(non_syn_ids, header=None, names=['refseqid', 'variant'], sep=' ')
queries = non_syn_vars.refseqid.unique()
records = parse_fasta(non_syn_fa)

with (wd_path / report_name).open('wt') as fout_report, \
        (wd_path / excluded_name).open('wt') as fout_excluded, \
        non_syn_fa_filtered.open('wt') as fout_seqs:
    failed, missing_seq, error_length, error_index, error_mismatch, error_reverse = {}, [], 0, 0, 0, 0
    for refseqid in queries:
        variants = non_syn_vars[non_syn_vars.refseqid == refseqid].variant
        seq = records.get(refseqid, None)
        if seq is None:
            missing_seq += [refseqid]
            continue
        seq_len = len(seq)
        violates_length_limit = seq_len > SNAP_LIMIT
        for idx, v in zip(variants.index, variants):
            reason = None
            wt, pos, mut = v[0], int(v[1:-1]), v[-1]
            if violates_length_limit:
                reason = 'length'
                error_length += 1
            elif pos > seq_len or wt != seq[pos - 1]:
                if pos > seq_len:
                    reason = 'index'
                    error_index += 1
                elif mut == seq[pos - 1]:
                    reason = 'reverse'
                    error_reverse += 1
                else:
                    reason = 'mismatch'
                    error_mismatch += 1
            if reason is not None:
                non_syn_vars.drop(idx, inplace=True)
                failed[refseqid] = failed.get(refseqid, []) + [(v, reason)]
    remaining_refseqids = set(non_syn_vars.refseqid)
    out = [
        ('# variants excluded from varidb query', '\n'),
        ('variants excluded   :', error_length + error_index + error_mismatch + error_reverse),
        ('affected sequences  :', len(failed)),
        ('removed sequences   :', len(records) - len(remaining_refseqids)),
        ('remaining sequences :', len(remaining_refseqids)),
        ('\n------ REASONS ------', '\n'),
        (' - length   :', error_length),
        (' - index    :', error_index),
        (' - mismatch :', error_mismatch),
        (' - reverse  :', error_reverse),
    ] + ([(' - missing  seq :', len(missing_seq))] if missing_seq else [])
    fout_report.writelines([f'{_[0]} {_[1]}\n' for _ in out])
    for seqid, variants in failed.items():
        for var, reason in variants:
            fout_excluded.write(f'{seqid}\t{var}\t{reason}\n')
    non_syn_vars.to_csv(non_syn_ids_filtered, index=False, header=False, sep=' ')
    for rid, seq in records.items():
        if rid in remaining_refseqids:
            fout_seqs.write(f'>{rid}\n{seq}\n')
