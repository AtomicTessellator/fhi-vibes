""" provides AttributeDict """

from collections import OrderedDict

# from hilde.helpers.warnings import warn


class MultiOrderedDict(OrderedDict):
    """A dict that can store multiple values for a key"""

    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super().__setitem__(key, value)


class AttributeDict(OrderedDict):
    """ Ordered dictionary with attribute access """

    def __getattr__(self, attr):
        if attr in self:
            return self[attr]
        # warn(f"Attribute {attr} not in dictionary, return None.", level=1)
        raise AttributeError(f"Attribute {attr} not in dictionary, return None.")

    def __dict__(self):
        return self.to_dict()

    def __repr__(self):
        return str(self.to_dict())

    def to_dict(self):
        """ (recursively) return plain python dictionary """
        rep = {}
        for key, val in self.items():
            if isinstance(val, AttributeDict):
                val = val.to_dict()
            rep.update({key: val})

        return rep

    def as_dict(self):
        """ return plain python dictionary (Fireworks compatibility) """
        return dict(self)
