from itertools import groupby


def list_group_by(ls, key_func):
    ls = sorted(ls, key=key_func)
    return [(k, list(g)) for k, g in groupby(ls, key_func)]
