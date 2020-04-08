"""helpers for aims"""

from vibes.filenames import filenames
from vibes.helpers import talk, warn


def peek_aims_uuid(file: str = filenames.output.aims) -> str:
    """peek into aims.out and find the uuid

    Args:
        file: Path to the aims.out file (default: aims.out)

    Returns:
        The uuid of an FHI-Aims calculation

    """
    try:
        with open(file) as f:
            line = next(l for l in f if "aims_uuid" in l)
            return line.split()[-1]
    except FileNotFoundError:
        talk(f"No aims_uuid found.")
    except StopIteration:
        warn(f"{file} presumably empty, no aims_uuid found")
    return ""


def get_aims_uuid_dict(file: str = filenames.output.aims) -> dict:
    """return aims uuid as dictionary

    Args:
        file: Path to the aims.out file (default: aims.out)

    Returns:
        The uuid of an FHI-Aims calculation as a dict

    """
    uuid = peek_aims_uuid(file)

    if uuid:
        return {"aims_uuid": uuid}
    return {}
