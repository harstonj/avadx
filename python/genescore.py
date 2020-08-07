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
        return self.gene.score_gene(series)

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


def process_snap_scores(row):
    res = np.nan
    score = row['snapfun.score']
    var_class = row['variant_class']
    is_na = np.isnan(score)
    if is_na and var_class == 'nonsynonymous SNV':
        res = np.nan
    elif is_na and var_class in ['stoploss', 'stopgain']:
        res = 100
    elif is_na and var_class == 'synonymous SNV':
        res = -99
    else:
        res = score
    return res


def extract_neutrality(row):
    label = np.nan
    score = row['snapfun.score']
    var_class = row['variant_class']
    if score > 0 and var_class == 'nonsynonymous SNV':
        label = 'NonNeutral'
    elif score <= 0 and var_class == 'nonsynonymous SNV':
        label = 'Neutral'
    else:
        label = np.nan
    return label


def run(annotations, scoretable, metadata, tools, variantfn, genefn, normalize, out):
    sample_name = Path(annotations.name.replace('.exonic_variant_function', '')).stem

    # load scoring functions
    scoring_functions = ScoringFunction(variantfn, genefn)

    # parse protein length info file
    protL = pd.read_csv(metadata)

    # parser annotations file
    exonic = pd.read_table(annotations, header=None)

    score_table = pd.DataFrame(columns=['Transcript', 'Mut'])
    for tool in tools:
        scoretable_file = (scoretable.parent / f'{scoretable.stem}_{tool}{scoretable.suffix}') if len(tools) > 1 else scoretable

        # parse score table file
        tool_score_table = pd.read_table(scoretable_file).rename(columns={'query': 'Transcript', 'variant': 'Mut'})
        score_table = pd.merge(score_table, tool_score_table, on=['Transcript', 'Mut'], how='outer')  # join with protein length information

    # preprocess annotations
    exonic_flted = exonic[[1, 2, 8]]  # subset columns of interest
    exonic_flted.columns = ['variant_class', 'variants_detailed', 'zygosity']  # rename columns
    exonic_flted = exonic_flted[exonic_flted.variants_detailed != 'UNKNOWN']  # drop all where annotations is not known
    exonic_flted = exonic_flted.assign(variants_detailed=exonic_flted.variants_detailed.str.split(',')).explode('variants_detailed')  # split multiple annotations into single rows
    exonic_flted = exonic_flted.dropna(subset=['variants_detailed'])  # drop all rows with missing variants_detailed caused by explode
    exonic_flted[['Gene', 'Transcript', 'exon', 'nuChange', 'aaChange']] = exonic_flted.variants_detailed.str.split(':', expand=True)  # seperate annotation details into rows
    exonic_flted = exonic_flted[exonic_flted.Gene != '']  # drop all rows with missing Gene caused by explode
    exonic_flted = exonic_flted.drop('variants_detailed', axis=1)  # drop old variants_detailed column
    exonic_flted = pd.merge(exonic_flted, protL[['Transcript', 'Prot_length']], on='Transcript', how='left')  # join with protein length information
    exonic_flted = exonic_flted.dropna(subset=['Prot_length'])  # drop all rows with missing Prot_length
    exonic_flted = exonic_flted.assign(Mut=exonic_flted.aaChange.str.replace('p.', ''))  # format variants
    exonic_flted = pd.merge(exonic_flted, score_table, on=['Transcript', 'Mut'], how='left')  # join with scores table

    # adjust scores
    exonic_flted['snapfun.score'] = exonic_flted.apply(lambda row: process_snap_scores(row), axis=1)  # process SNAPfun scores
    exonic_flted = exonic_flted.dropna(subset=['snapfun.score'])  # drop all rows with missing snapfun.score

    # Generate a summary for this sample:
    exonic_flted_summary = exonic_flted[['variant_class', 'zygosity', 'snapfun.score']].copy()
    exonic_flted_summary['neutrality'] = exonic_flted_summary.apply(lambda row: extract_neutrality(row), axis=1)
    summary_exonic = exonic_flted_summary.groupby(['variant_class', 'zygosity', 'neutrality']).count().rename(columns={'snapfun.score': 'counts'}).reset_index()
    summary_exonic.to_csv(out / f'summary_{sample_name}.csv', sep=',', index=False)

    # NOTE def generate_genescores():

    df = exonic_flted.copy()

    # compute variant scores
    df['variant_score'] = df.apply(scoring_functions.score_variant, axis=1)

    # aggregate single variant scores to per gene score (gene_score)
    df_gs = df.groupby(['Gene', 'Transcript']).agg({'variant_score': scoring_functions.score_gene, 'Prot_length': 'first'}).rename(columns={'variant_score': 'gene_score'}).reset_index()

    # normalize variant score by protein lengths
    df_gs['gene_score_normalized'] = df_gs[['gene_score', 'Prot_length']].apply(lambda x: x['gene_score'] * 100 / x['Prot_length'], axis=1)

    # set precision
    df_gs = df_gs.round({'gene_score': 4, 'gene_score_normalized': 4})

    if normalize == 'n':
        df_gs.to_csv(out / f'{sample_name}.gs', columns=['Gene', 'Transcript', 'gene_score'], index=False)
    elif normalize == 'y':
        df_gs.to_csv(out / f'{sample_name}.gs', columns=['Gene', 'Transcript', 'gene_score_normalized'], index=False)
    elif normalize == 'both':
        df_gs.to_csv(out / f'{sample_name}.gs', columns=['Gene', 'Transcript', 'gene_score', 'gene_score_normalized'], index=False)


def main(args):
    run(args.annotations, args.scoretable, args.metadata, [_.strip() for _ in args.tools.split(',')], args.variantfn, args.genefn, args.normalize, args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--annotations', type=Path,
        help='Single sample ANNOVAR annotations (.exonic_variant_function)'
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
        '-n', '--normalize', type=str, default='both', choices=['y', 'n', 'both'],
        help='compute absolute [n], normalized [y] or both [both] gene scores. Normalization is [0-1] using: (gene_score * 100 / protein_length)'
    )
    parser.add_argument(
        '-o', '--out', type=Path,
        help='path to the output folder'
    )
    main(parser.parse_args())
