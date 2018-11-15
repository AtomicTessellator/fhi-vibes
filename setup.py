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
    author_email="knoop@fhi-berlin.mpg.de, purcell@fhi-berlin.mpg.de",
    license="Unlicense",
    # install_requires=['numpy', 'phonopy', 'spglib'],
    packages=setuptools.find_packages(),
    ext_modules=[ext],
    scripts=[
        "hilde/scripts/geometry_info",
        "hilde/scripts/get_relaxation_info",
        "hilde/scripts/make_supercell",
        "hilde/scripts/refine_geometry",
    ],
    entry_points={
        'console_scripts': [
            'qlaunch_hilde = hilde.fireworks_api_adapter.scripts.qlaunch_run:qlaunch',
            'rlaunch_hilde = hilde.fireworks_api_adapter.scripts.rlaunch_run:rlaunch',
            'claunch_hilde = hilde.fireworks_api_adapter.scripts.claunch_run:claunch'
        ]
    },
    zip_safe=False,
)
