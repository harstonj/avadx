import os
import argparse
import importlib
import pandas as pd
import numpy as np
from pathlib import Path


class ScoringFunction:
    def __init__(self, variant='default_variant', gene='default_gene_sum'):
        self.import_modules(name='scoring_functions', variant=Path(variant), gene=Path(gene))

    def score_variant(self, row):
        return self.variant.score_variant(row)

    def score_gene(self, series):
        series_missing_removed = series.dropna()
        aggregated_score = np.nan if series_missing_removed.empty else self.gene.score_gene(series_missing_removed)
        return aggregated_score

    def import_modules(self, name, package: str = None, variant: Path = None, gene: Path = None):
        variant_file = variant if variant else Path(variant.replace('.', os.sep))
        gene_file = gene if gene else Path(gene.replace('.', os.sep))
        if not variant_file.exists():
            variant_file = Path(__file__).parent.absolute() / name / f'{variant}.py'
        if not gene_file.exists():
            gene_file = Path(__file__).parent.absolute() / name / f'{gene}.py'
        variant_spec = importlib.util.spec_from_file_location(f'scoring_functions.{variant}', variant_file)
        self.variant = importlib.util.module_from_spec(variant_spec)
        variant_spec.loader.exec_module(self.variant)
        gene_spec = importlib.util.spec_from_file_location(f'scoring_functions.{gene}', gene_file)
        self.gene = importlib.util.module_from_spec(gene_spec)
        gene_spec.loader.exec_module(self.gene)


def extract_neutrality(row):
    label = np.nan
    score = row['snapfun.score']
    var_type = row['type']
    if score > 0 and var_type == 'nonsynonymous SNV':
        label = 'NonNeutral'
    elif score <= 0 and var_type == 'nonsynonymous SNV':
        label = 'Neutral'
    else:
        label = np.nan
    return label


def run(annotations, indels, scoretable, metadata, tools, variantfn, genefn, normalize, out):
    sample_name = Path(annotations.name.replace('.exonic_variant_function', '')).stem

    # load scoring functions
    scoring_functions = ScoringFunction(variantfn, genefn)

    # parse protein length info file
    protL = pd.read_csv(metadata).rename(columns={'Transcript': 'transcript', 'Prot_length': 'prot_length'})

    # parser annotations file
    exonic_snp = pd.read_table(annotations, header=None)
    exonic_snp['class'] = 'snp'
    if indels:
        indel_annotations = annotations.parent.parent / (annotations.parent.name + '_indels') / annotations.name
        exonic_indels = pd.read_table(indel_annotations, header=None)
        exonic_indels['class'] = 'indel'
        exonic_all = pd.concat([exonic_snp, exonic_indels], ignore_index=True)
    else:
        exonic_all = exonic_snp

    score_table = pd.DataFrame(columns=['transcript', 'variant'])
    for tool in tools:
        scoretable_file = (scoretable.parent / f'{scoretable.stem}_{tool}{scoretable.suffix}') if len(tools) > 1 else scoretable

        # parse score table file
        tool_score_table = pd.read_table(scoretable_file).rename(columns={'query': 'transcript', 'variant': 'variant'})
        score_table = pd.merge(score_table, tool_score_table, on=['transcript', 'variant'], how='outer')  # join with protein length information

    # preprocess annotations
    exonic_flted = exonic_all[[1, 2, 8, 'class']]  # subset columns of interest
    exonic_flted.columns = ['type', 'detail', 'zygosity', 'class']  # rename columns
    # NOTE: removed filtering of unknowns - those are set to NaN or an actual value using the plugin scoring function
    # exonic_flted = exonic_flted[exonic_flted.detail != 'UNKNOWN']  # drop all where annotations is not known
    exonic_flted = exonic_flted.assign(detail=exonic_flted.detail.str.split(',')).explode('detail')  # split multiple annotations into single rows
    exonic_flted = exonic_flted.dropna(subset=['detail'])  # drop all rows with missing detail caused by explode
    exonic_flted[['gene', 'transcript', 'exon', 'nu_change', 'aa_change']] = exonic_flted.detail.str.split(':', expand=True)  # seperate annotation details into rows
    exonic_flted = exonic_flted[exonic_flted.gene != '']  # drop all rows with missing gene caused by explode
    exonic_flted = exonic_flted.drop('detail', axis=1)  # drop old detail column
    exonic_flted = pd.merge(exonic_flted, protL[['transcript', 'prot_length']], on='transcript', how='left')  # join with protein length information
    exonic_flted = exonic_flted.dropna(subset=['prot_length'])  # drop all rows with missing prot_length
    exonic_flted = exonic_flted.assign(variant=exonic_flted.aa_change.str.replace('p.', ''))  # format variants
    exonic_flted = pd.merge(exonic_flted, score_table, on=['transcript', 'variant'], how='left')  # join with scores table

    # Generate a summary for this sample:
    exonic_flted_summary = exonic_flted[['type', 'zygosity', 'snapfun.score']].copy()
    exonic_flted_summary['neutrality'] = exonic_flted_summary.apply(lambda row: extract_neutrality(row), axis=1)
    summary_exonic = exonic_flted_summary.groupby(['type', 'zygosity', 'neutrality']).count().rename(columns={'snapfun.score': 'counts'}).reset_index()
    summary_exonic.to_csv(out / f'summary_{sample_name}.csv', sep=',', index=False)

    # copy score table
    df = exonic_flted.copy()

    # compute variant scores
    df['variant_score'] = df.apply(scoring_functions.score_variant, axis=1)

    # aggregate single variant scores to per gene score (gene_score)
    df_gs = df.groupby(['gene', 'transcript']).agg({'variant_score': scoring_functions.score_gene, 'prot_length': 'first'}).rename(columns={'variant_score': 'gene_score'}).reset_index()

    # normalize variant score by protein lengths
    df_gs['gene_score_normalized'] = df_gs[['gene_score', 'prot_length']].apply(lambda x: x['gene_score'] * 100 / x['prot_length'], axis=1)

    # set precision and rename columns
    df_gs = df_gs.round({'gene_score': 4, 'gene_score_normalized': 4}).rename(columns={'gene': 'Gene', 'transcript': 'Transcript'})

    # save genescores
    df_gs.to_csv(out / f'{sample_name}.gs', columns=['Gene', 'Transcript', 'gene_score', 'gene_score_normalized'], index=False)


def merge(results, metadata, cvscheme, variantfn, genefn, normalize, out):
    # load cv-scheme
    cvscheme_df = pd.read_csv(cvscheme, header=None, usecols=[0], names=['sampleid'], dtype={'sampleid': str})
    cvscheme_df = cvscheme_df.set_index('sampleid')

    # merge sample genescore files
    genescores_merged = pd.DataFrame(columns=['Gene', 'Transcript'])
    genescores_samples = {_.name.replace('sample.', '').replace('.gs', ''): _ for _ in results.glob('*.gs')}
    for sample_id, genescores_file in genescores_samples.items():
        genescores_merged = genescores_merged.merge(pd.read_csv(genescores_file, dtype={'gene_score': np.float64, 'gene_score_normalized': np.float64}).rename(columns={'gene_score': f'{sample_id}_r', 'gene_score_normalized': f'{sample_id}_n'}), on=['Gene', 'Transcript'], how='outer')

    # load scoring functions
    scoring_functions = ScoringFunction(variantfn, genefn)

    # parse protein length info file
    protL = pd.read_csv(metadata).rename(columns={'Transcript': 'transcript', 'Prot_length': 'prot_length'})

    # filter by and save raw and/pr normalized scores and process missing values
    if normalize == 'no':
        merged_out_path = out / 'GeneScoreTable_raw.csv'
        genescores_merged_r = genescores_merged.filter(regex='Gene|Transcript|_r$', axis=1)
        genescores_merged_r.columns = genescores_merged_r.columns.str.rstrip('_r')
        genescores_merged_r = genescores_merged_r.fillna(scoring_functions.gene.NA_SCORE)
        genescores_merged_r[['Gene', 'Transcript'] + cvscheme_df.index.tolist()].to_csv(merged_out_path, index=False, float_format='%.4g')
    elif normalize == 'yes':
        merged_out_path = out / 'GeneScoreTable_normalized.csv'
        genescores_merged_n = genescores_merged.filter(regex='Gene|Transcript|_n$', axis=1)
        genescores_merged_n.columns = genescores_merged_n.columns.str.rstrip('_n')
        genescores_merged_n = pd.merge(genescores_merged_n, protL[['transcript', 'prot_length']], left_on='Transcript', right_on='transcript', how='left').drop('transcript', axis=1)
        genescores_merged_n = genescores_merged_n.apply(lambda x: x.fillna(scoring_functions.gene.NA_SCORE / x['prot_length']), axis=1).drop('prot_length', axis=1)
        genescores_merged_n[['Gene', 'Transcript'] + cvscheme_df.index.tolist()].to_csv(merged_out_path, index=False, float_format='%.4g')


def main(args):
    if args.merge:
        merge(args.results, args.metadata, args.cvscheme, args.variantfn, args.genefn, args.normalize, args.out)
    else:
        run(args.annotations, args.indels, args.scoretable, args.metadata, [_.strip() for _ in args.tools.split(',')], args.variantfn, args.genefn, args.normalize, args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-M', '--merge', action='store_true',
        help='merge single sample genescore (.gs) files'
    )
    parser.add_argument(
        '-r', '--results', type=Path,
        help='path to folder containing all single sample genescores (.gs) files'
    )
    parser.add_argument(
        '-a', '--annotations', type=Path,
        help='single sample ANNOVAR annotations (.exonic_variant_function)'
    )
    parser.add_argument(
        '-i', '--indels', action='store_true',
        help='also use ANNOVAR annotations for indels (.exonic_variant_function)'
    )
    parser.add_argument(
        '-s', '--scoretable', type=Path,
        help='varidb score table'
    )
    parser.add_argument(
        '-m', '--metadata', type=Path,
        help='metadata: protein length information (RefSeq ids)'
    )
    parser.add_argument(
        '-c', '--cvscheme', type=Path,
        help='path to cross-validation scheme file'
    )
    parser.add_argument(
        '-t', '--tools', type=str, default='SNAPfun',
        help='comma separated list of prediction methods to use for scoring'
    )
    parser.add_argument(
        '-g', '--genefn', type=str, default='default_gene_sum',
        help='module name to use for aggregating single variant scores into final gene score'
    )
    parser.add_argument(
        '-v', '--variantfn', type=str, default='default_variant',
        help='module name to compute single variant score'
    )
    parser.add_argument(
        '-n', '--normalize', type=str, default='no', choices=['yes', 'no'],
        help='compute raw [no] or normalized [y] gene scores. Normalization is [0-1] using: (gene_score * 100 / protein_length)'
    )
    parser.add_argument(
        '-o', '--out', type=Path,
        help='path to the output folder'
    )
    main(parser.parse_args())
