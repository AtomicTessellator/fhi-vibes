"""helpers for lists and tuples"""

from itertools import groupby


def reduce_list(obj, reduce=True):
    """Reduce a with duplicate entries and return tuples of (count, entry)"""
    if reduce:
        return tuple((len(list(g)), k) for k, g in groupby(obj))
    return obj


def expand_list(obj):
    """Expand a list of tuples (count, entry) as produced ty `reduce_list`"""
    if isinstance(obj[0], type(obj)):
        lis = []
        for l in (int(l) * [g] for (l, g) in obj):
            lis.extend(l)
        return lis
    return obj


def list_dim(a: list) -> int:
    """Dimension of a (nested) pure python list, similar to np.shape"""
    if not isinstance(a, list):
        return []
    if a == []:
        return 0
    return [len(a), *list_dim(a[0])]


def list2ndr(lis: list) -> str:
    """Convert list to string"""
    return "[{}]".format(", ".join([str(el) for el in lis]))
