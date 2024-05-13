"""provide numpy extension to compile fortran routines"""

from numpy.distutils.core import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="vibes.helpers.supercell.supercell",
            sources=[
                "vibes/helpers/supercell/linalg.f90",
                "vibes/helpers/supercell/supercell.f90",
            ],
            extra_compile_args=["-O3"],
        )
    ]
)
