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
    version="0.1.3",
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
            "hilde_phonopy = hilde.scripts.hilde_phonopy:main",
            "create_samples= hilde.scripts.create_samples:main",
            "hilde_single_task= hilde.scripts.hilde_single_task:main",
            "hilde_run = hilde.scripts.hilde_run:main",
            "md_sum= hilde.scripts.md_sum:main",
            "yaml2json = hilde.scripts.yaml2json:main",
            "refine_geometry = hilde.scripts.refine_geometry:main",
            "add_workflow_to_lp = hilde.scripts.add_fireworks_workflow:main",
            "add_phonon_workflow_to_lp = hilde.scripts.add_fireworks_phonon_workflow:main",
            "add_traj_to_db = hilde.scripts.add_trajectory_to_database:main",
            "qlaunch_hilde = hilde.fireworks.scripts.qlaunch_run:qlaunch",
            "rlaunch_hilde = hilde.fireworks.scripts.rlaunch_run:rlaunch",
            "claunch_hilde = hilde.fireworks.scripts.claunch_run:claunch",
            "nomad_upload = hilde.scripts.nomad_upload:main",
            "update_md_trajectory = hilde.scripts.update_md_trajectory:main",
            "trajectory2tdep = hilde.scripts.trajectory2tdep:main",
            "trajectory2xyz = hilde.scripts.trajectory2xyz:main",
            "suggest_k_grid = hilde.scripts.suggest_k_grid:main",
        ]
    },
    data_files=[(".", ["hilde.cfg", "fireworks.cfg"])],
    zip_safe=False,
    install_requires=["numpy", "scipy", "ase"],
)
