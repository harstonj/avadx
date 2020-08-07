from functools import reduce


def score_gene(series):
    # score correlates positive with function change
    # score = 0: protein function disrupted (e.g. stopgain/loss variant)
    # score = 1: protein function unchanged (e.g. synonymous / nonsynonymous (neutral) variant

    return reduce(lambda x, y: x * y, 1 - series)
