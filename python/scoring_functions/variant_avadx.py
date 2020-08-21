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
                        'type', 'zygosity', 'class', 'gene', 'transcript', 'exon', 'nu_change',
                        'aa_change', 'prot_length', 'variant'
                    ]

                    Details for multi-value keys:
                     - type:
                        [
                            frameshift insertion,
                            frameshift deletion,
                            frameshift block substitution,
                            stopgain,
                            stoploss,
                            nonframeshift insertion,
                            nonframeshift deletion,
                            nonframeshift block substitution,
                            nonsynonymous SNV,
                            synonymous SNV,
                            unknown
                        ]
                     - zygosity:
                        [ het, hom ]
                     - class:
                        [ snp, indel ]

    Returns
    -------
    variant_score : float or None

                    Returned variant_score has to be standardized to range [0, 1]
                    If variant_score is None it will be ignored from further consideration in gene_score aggregation step
    """

    variant_score = None
    var_type = var_meta_dict['type']
    var_class = var_meta_dict['class']
    snapfun_score = var_meta_dict['snapfun.score']
    heti = .25 if var_meta_dict['zygosity'] == 'het' else 1

    # scoring for "snp" class variants
    if var_class == 'snp':
        if var_type == 'stopgain':
            variant_score = 1 * heti
        elif var_type == 'stoploss':
            variant_score = 1 * heti
        elif var_type == 'nonsynonymous SNV':
            if snapfun_score <= 0:
                variant_score = .055 * heti
            elif snapfun_score > 0:
                variant_score = (.06 + (snapfun_score / 100 * .94)) * heti
            else:
                variant_score = .055 * heti
        elif var_type == 'synonymous SNV':
            variant_score = .05 * heti

    # scoring for "indel" class variants
    elif var_class == 'indel':
        if var_type == 'frameshift insertion':
            variant_score = None
        elif var_type == 'frameshift deletion':
            variant_score = None
        elif var_type == 'frameshift block substitution':
            variant_score = None
        elif var_type == 'stopgain':
            variant_score = None
        elif var_type == 'stoploss':
            variant_score = None
        elif var_type == 'nonframeshift insertion':
            variant_score = None
        elif var_type == 'nonframeshift deletion':
            variant_score = None
        elif var_type == 'nonframeshift block substitution':
            variant_score = None
        elif var_type == 'nonsynonymous SNV':
            variant_score = None
        elif var_type == 'synonymous SNV':
            variant_score = None
        elif var_type == 'unknown':
            variant_score = None

    return variant_score
