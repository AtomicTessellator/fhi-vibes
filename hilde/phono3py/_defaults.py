""" hilde Phono3py defaults """

from hilde.phonopy import defaults as defaults_fc2

defaults = defaults_fc2.copy()
defaults.update({"is_diagonal": True, "cutoff_pair_distance": 100.0, "log_level": 2})
