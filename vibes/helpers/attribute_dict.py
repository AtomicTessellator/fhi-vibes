""" provides AttributeDict """

from collections import OrderedDict


class AttributeDict(OrderedDict):
    """ Ordered dictionary with attribute access """

    def __getattr__(self, attr):
        if attr in self:
            return self[attr]
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
