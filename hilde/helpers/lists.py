"""helpers for lists"""


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

