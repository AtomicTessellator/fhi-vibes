"""helpers for lists and tuples"""
from itertools import groupby


def reduce_list(obj):
    """reduce a with duplicate entries and return tuples of (count, entry)"""
    return tuple((len(list(g)), k) for k, g in groupby(obj))


def expand_list(obj):
    """expand a list of tuples (count, entry) as produced ty `reduce_list`"""
    if isinstance(obj[0], type(obj)):
        lis = []
        for l in (l * [g] for (l, g) in obj):
            lis.extend(l)
        return lis
    return obj


def list_dim(a):
    """dimension of a (nested) pure Python list

    Parameters
    ----------
    a: list
        The input list

    Returns
    -------
    int
        The dimension of the pure python list
    """
    if not type(a) == list:
        return []
    return [len(a)] + list_dim(a[0])


def list2str(lis):
    """convert list to string

    Parameters
    ----------
    lis: list
        list to convert to str

    Returns
    -------
    str
        The json string version of the list
    """
    return "[{}]".format(", ".join([str(el) for el in lis]))

