import collections

# Avalanche
Fmax = 100

# file formats
_dct = {"geometry": "aims"}
suffixes = collections.namedtuple("format", _dct.keys())(**_dct)
