import setuptools
from numpy.distutils.core import setup, Extension

ext = Extension(
    name="hilde.helpers.supercell.supercell",
    sources=[
        "hilde/helpers/supercell/linalg.f90",
        "hilde/helpers/supercell/supercell.f90",
    ],
    extra_compile_args=["-O3"],
)

setup(
    name="HiLDe",
    version="0.1",
    description="Haber Institute Lattice Dynamics Package",
    author="Florian Knoop, Thomas Purcell",
    author_email="knoop@fhi-berlin.mpg.de",
    license="MIT License",
    # install_requires=['numpy', 'phonopy', 'spglib'],
    packages=setuptools.find_packages(),
    ext_modules=[ext],
    scripts=[
        "hilde/scripts/geometry_info",
        "hilde/scripts/get_relaxation_info",
        "hilde/scripts/make_supercell",
        "hilde/scripts/refine_geometry",
    ],
    zip_safe=False,
)
