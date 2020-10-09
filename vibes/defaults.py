import collections


# Avalanche
F_max = 10
F_window = 50

# file formats
_dct = {"geometry": "aims"}
format = collections.namedtuple("format", _dct.keys())(**_dct)
