"""provide numpy extension to compile fortran routines"""
from numpy.distutils.core import setup, Extension  # noqa: F401

ext = Extension(
    name="vibes.helpers.supercell.supercell",
    sources=[
        "vibes/helpers/supercell/linalg.f90",
        "vibes/helpers/supercell/supercell.f90",
    ],
    extra_compile_args=["-O3"],
)


def build(setup_kwargs):
    """add Extension to setup kwargs"""

    setup_kwargs.update({"ext_modules": [ext]})
