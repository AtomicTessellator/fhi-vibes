""" provides AttributeDict """

from collections import OrderedDict
from hilde.helpers.warnings import warn


class AttributeDict(OrderedDict):
    def __getattr__(self, attr):
        if attr in self:
            return self[attr]
        else:
            warn(f"Attribute {attr} not in dictionary, return None.", level=1)
            return None

    def __setattr__(self, attr, value):
        self[attr] = value

    def to_dict(self):
        """ return plain python dictionary """
        return dict(self)

    def as_dict(self):
        """ return plain python dictionary (Fireworks compatibility) """
        return dict(self)
