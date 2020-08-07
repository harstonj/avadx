from functools import reduce


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

    return reduce(lambda x, y: x + y, variant_scores)
