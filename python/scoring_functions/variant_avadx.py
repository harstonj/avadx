import numpy as np


def score_variant(row):
    zygosity_het = .25
    zygosity_hom = 1

    var_class = row['variant_class']
    zygosity_class = row['zygosity']
    snapfun_score = row['snapfun.score']
    heti = zygosity_het if zygosity_class == 'het' else zygosity_hom
    variant_score = np.nan

    if np.isnan(snapfun_score):
        if var_class == 'nonsynonymous SNV':
            snapfun_score = np.nan
        elif var_class in ['stoploss', 'stopgain']:
            snapfun_score = 100
        elif var_class == 'synonymous SNV':
            snapfun_score = -99

    if var_class == 'synonymous SNV':
        variant_score = 0.05 * heti
    elif var_class == 'nonsynonymous SNV' and snapfun_score <= 0:
        variant_score = 0.055 * heti
    elif snapfun_score > 0:
        variant_score = (0.06 + (snapfun_score / 100 * .94)) * heti
    else:
        variant_score = 0.055 * heti

    return variant_score
