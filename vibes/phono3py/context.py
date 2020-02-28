"""Phono3py workflow context managing"""

from vibes.phonopy.context import PhonopyContext

from . import postprocess, wrapper
from ._defaults import defaults, mandatory, name


class Phono3pyContext(PhonopyContext):
    """PhonopyContext with changed name"""

    def __init__(self, *args, **kwargs):
        kwargs.update(
            {"name": name, "defaults_kw": defaults, "mandatory_kw": mandatory}
        )
        super().__init__(*args, **kwargs)

        self.backend = wrapper
        self.postprocess = postprocess
