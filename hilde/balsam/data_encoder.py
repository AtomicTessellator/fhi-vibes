"""JSONable encoder/decoder for the data object Modified from ase.io.jsonio"""
import datetime
import json

import numpy as np


class MyEncoder(json.JSONEncoder):
    """Encoder taken from ASE to encode numpy objects into JSONable objects"""

    def default(self, obj):
        """default encode function for obj"""
        if hasattr(obj, "todict") or hasattr(obj, "to_dict"):
            try:
                d = obj.todict()
            except AttributeError:
                d = obj.to_dict()

            if not isinstance(d, dict):
                raise RuntimeError(
                    "todict() of {} returned object of type {} "
                    "but should have returned dict".format(obj, type(d))
                )

            return d
        if isinstance(obj, np.ndarray) or hasattr(obj, "__array__"):
            if obj.dtype == complex:
                return {"__complex_ndarray__": (obj.real.tolist(), obj.imag.tolist())}
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, datetime.datetime):
            return {"__datetime__": obj.isoformat()}
        return json.JSONEncoder.default(self, obj)


encode = MyEncoder().encode


def object_hook(dct):
    """Get an object from the JSON encoded dictionary"""
    if "__datetime__" in dct:
        return datetime.datetime.strptime(dct["__datetime__"], "%Y-%m-%dT%H:%M:%S.%f")
    if "__complex_ndarray__" in dct:
        r, i = (np.array(x) for x in dct["__complex_ndarray__"])
        return r + i * 1j

    if "__ase_objtype__" in dct:
        objtype = dct.pop("__ase_objtype__")
        dct = numpyfy(dct)

        # We just try each object type one after another and instantiate
        # them manually, depending on which kind it is.
        # We can formalize this later if it ever becomes necessary.
        if objtype == "cell":
            from ase.geometry.cell import Cell

            obj = Cell(**dct)
        elif objtype == "bandstructure":
            from ase.dft.band_structure import BandStructure

            obj = BandStructure(**dct)
        elif objtype == "bandpath":
            from ase.dft.kpoints import BandPath

            obj = BandPath(**dct)
        else:
            raise RuntimeError(
                "Do not know how to decode object type {} "
                "into an actual object".format(objtype)
            )

        assert obj.ase_objtype == objtype
        return obj

    return dct


mydecode = json.JSONDecoder(object_hook=object_hook).decode


def intkey(key):
    """get the integer of the key"""
    try:
        return int(key)
    except ValueError:
        return key


def numpyfy(obj):
    """Create a numpy object from the JSONAble object"""
    if isinstance(obj, dict):
        if "__complex_ndarray__" in obj:
            r, i = (np.array(x) for x in obj["__complex_ndarray__"])
            return r + i * 1j
        return dict((intkey(key), numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            pass
        else:
            if a.dtype in [bool, int, float]:
                return a
        obj = [numpyfy(value) for value in obj]
    return obj


def decode(txt):
    """Decode the JSON string into a dictionary"""
    return numpyfy(mydecode(txt))
