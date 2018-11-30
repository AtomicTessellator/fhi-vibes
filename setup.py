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
    entry_points={
        "console_scripts": [
            "geometry_info = hilde.scripts.geometry_info:main",
            "rewrite_geometry= hilde.scripts.rewrite_geometry:main",
            "get_relaxation_info = hilde.scripts.get_relaxation_info:main",
            "make_supercell = hilde.scripts.make_supercell:main",
            "phonopy_info= hilde.scripts.phonopy_info:main",
            "create_samples= hilde.scripts.create_samples:main",
            "hilde_single_task= hilde.scripts.hilde_single_task:main",
            "md_sum= hilde.scripts.md_sum:main",
            "yaml2json= hilde.scripts.yaml2json:main",
            "refine_geometry = hilde.scripts.refine_geometry:main",
            "add_workflow_to_lp = hilde.scripts.add_fireworks_workflow:main",
            "qlaunch_hilde = hilde.fireworks_api_adapter.scripts.qlaunch_run:qlaunch",
            "rlaunch_hilde = hilde.fireworks_api_adapter.scripts.rlaunch_run:rlaunch",
            "claunch_hilde = hilde.fireworks_api_adapter.scripts.claunch_run:claunch",
            "nomad_upload = hilde.scripts.nomad_upload:main",
        ]
    },
    data_files=[(".", ["hilde.cfg"])],
    zip_safe=False,
)
