#!/usr/bin/env python

import setuptools
from numpy.distutils.core import setup, Extension

ext = Extension(name = 'hilde.helpers.supercell.supercell',
                sources=['hilde/helpers/supercell/supercell.f90',
                         'hilde/helpers/supercell/linalg.f90'],
                extra_compile_args = ['-O3'],
)

setup(name='HiLDe',
      version='0.1',
      description='Haber Institute Lattice Dynamics Package',
      author='Thomas Purcell, Florian Knoop',
      author_email='knoop@fhi-berlin.mpg.de',
      license='MIT License',
      # install_requires=['numpy', 'phonopy', 'spglib'],
      packages=setuptools.find_packages(),
      ext_modules=[ext],
      scripts=['hilde/scripts/geometry_info']
      )
