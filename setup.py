#!/usr/bin/env python

from numpy.distutils.core import setup
from numpy.distutils.core import Extension

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
      install_requires=['numpy', 'phonopy', 'spglib'],
      packages=['hilde.helpers', 'hilde.konstanten', 'hilde.parsers',
                'hilde.phonopy', 'hilde.settings', 'hilde.structure',
                'hilde.tasks', 'hilde.helpers.supercell'],
      ext_modules=[ext]
      )
