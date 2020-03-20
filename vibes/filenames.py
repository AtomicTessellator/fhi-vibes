"""default filenames"""
import collections
from pathlib import Path

from vibes import keys
from vibes._defaults import DEFAULT_GEOMETRY_FILE

_filename = Path(__file__).stem

_geometry = DEFAULT_GEOMETRY_FILE

_suffixes = "suffixes"
_dct = {"next_step": ".next_step"}
suffixes = collections.namedtuple(_suffixes, _dct.keys())(**_dct)


_output = "output"
_dct = {"aims": "aims.out"}
output = collections.namedtuple(_output, _dct.keys())(**_dct)


_dct = {
    "atoms": _geometry,
    "atoms_next": _geometry + suffixes.next_step,
    "primitive": _geometry + ".primitive",
    "supercell": _geometry + ".supercell",
    "trajectory": keys.trajectory + ".son",
    "trajectory_dataset": keys.trajectory + ".nc",
    _output: output,
}
filenames = collections.namedtuple(_filename, _dct.keys())(**_dct)
