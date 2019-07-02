"""helpers for aims"""

from hilde.helpers import talk


def peek_aims_uuid(filename="aims.out"):
    """peek into aims.out and find the uuid

    Parameters
    ----------
    filename: str
        Path to the aims.out file (default: aims.out)

    Returns
    -------
    str
        The uuid of an FHI-Aims calculation
    """
    try:
        with open(filename) as f:
            line = next(l for l in f if "aims_uuid" in l)
            return line.split()[-1]
    except FileNotFoundError:
        talk(f"No aims_uuid found.")
        return ""


def get_aims_uuid_dict(filename="aims.out"):
    """return aims uuid as dictionary

    Parameters
    ----------
    filename: str
        Path to the aims.out file (default: aims.out)

    Returns
    -------
    dict
        The uuid of an FHI-Aims calculation as a dict
    """
    uuid = peek_aims_uuid(filename)

    if uuid:
        return {"aims_uuid": uuid}
    return {}
