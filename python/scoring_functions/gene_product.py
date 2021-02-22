from functools import reduce
import numpy as np


# gene score if no variants scores are available
NA_SCORE = 1


def score_gene(variants_df):
    """ Alternative AVA,Dx variant score aggregation function (product)

    Parameters
    ----------
    variants_df :   pandas.DataFrame
                    DataFrame containing all variant scores and metadata for the current gene

                    metadata keys are:
                    [
                        'type', 'zygosity', 'class', 'gene', 'transcript', 'exon', 'nu_change',
                        'aa_change', 'prot_length', 'variant', 'snapfun.score', 'variant_score'
                    ]

    Returns
    -------
    gene_score : float
    """

    # get list of variant scores
    variant_scores = variants_df.variant_score

    # drop all missing scores
    variant_scores_valid = variant_scores.dropna()

    # return NA in case all scores are missing, otherwise (product)
    gene_score = np.nan if variant_scores_valid.empty else reduce(lambda x, y: x * y, 1 - variant_scores)

    return gene_score
