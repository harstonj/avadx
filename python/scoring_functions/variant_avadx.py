import numpy as np


def score_variant(row):
    zygosity_het = .25
    zygosity_hom = 1

    var_class = row['variant_class']
    zygosity_class = row['zygosity']
    snapfun = row['snapfun.score']
    heti = zygosity_het if zygosity_class == 'het' else zygosity_hom

    score = np.nan

    if var_class == 'synonymous SNV':
        score = 0.05 * heti
    elif var_class == 'nonsynonymous SNV' and snapfun <= 0:
        score = 0.055 * heti
    elif snapfun > 0:
        score = (0.06 + (snapfun / 100 * .94)) * heti
    else:
        score = np.nan

    return score
