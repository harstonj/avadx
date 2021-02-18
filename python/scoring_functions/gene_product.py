from functools import reduce


# gene score if no variants scores are available
NA_SCORE = 1


def score_gene(variant_scores, gene, transcript):
    """ Alternative AVA,Dx variant score aggregation function (product)

    Parameters
    ----------
    variant_scores : pandas.core.series.Series
                     Series contains all variant scores computed for the current gene
    gene           : gene name
    transcript     : transcript name

    Returns
    -------
    gene_score : float
    """

    return reduce(lambda x, y: x * y, 1 - variant_scores)
