# AVAD,X scoring function

def score_variant(var_meta_dict):
    """ Default AVA,Dx scoring function using SNAPfun scores

    Parameters
    ----------
    var_meta_dict : pandas.core.series.Series

                    Series contains all scores extracted from varidb (for SNAPfun and other selected tools):
                    [
                         'snapfun.score'
                    ]

                    Additional available keys are:
                    [
                        'type', 'zygosity', 'gene', 'transcript', 'exon', 'nu_change',
                        'aa_change', 'prot_length', 'variant'
                    ]

    Returns
    -------
    variant_score : float

                    Returned variant_score has to be standardized to range [0, 1]
    """

    var_type = var_meta_dict['type']
    snapfun_score = var_meta_dict['snapfun.score']
    heti = .25 if var_meta_dict['zygosity'] == 'het' else 1

    if var_type == 'synonymous SNV':
        variant_score = .05 * heti
    elif var_type == 'nonsynonymous SNV' and snapfun_score <= 0:
        variant_score = .055 * heti
    elif var_type in ['stoploss', 'stopgain']:
        variant_score = 1 * heti
    elif snapfun_score > 0:
        variant_score = (.06 + (snapfun_score / 100 * .94)) * heti
    else:
        variant_score = .055 * heti

    return variant_score
