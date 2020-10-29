from functools import reduce


# gene score if no variants scores are available
NA_SCORE = 0
MAX_SCORE = 3


def score_gene(variant_scores):
    """ Default AVA,Dx variant score aggregation function (sum)

    Parameters
    ----------
    variant_scores : pandas.core.series.Series
                     Series contains all variant scores computed for the current gene

    Returns
    -------
    gene_score : float
    """

    return max(reduce(lambda x, y: x + y, variant_scores), MAX_SCORE)
