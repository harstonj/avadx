from functools import reduce


def score_gene(variant_scores):
    """ Alternative AVA,Dx variant score aggregation function (product)

    Parameters
    ----------
    variant_scores : pandas.core.series.Series
                     Series contains all variant scores computed for the current gene

    Returns
    -------
    gene_score : float
    """

    return reduce(lambda x, y: x * y, 1 - variant_scores)
