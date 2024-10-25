project('fhi-vibes',
  'c', 'fortran',
  version: '1.1.0',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

supercell_source = custom_target('supercellmodule',
  input : [
    'vibes/helpers/supercell/supercell.f90',
    'vibes/helpers/supercell/linalg.f90',
    ],
  output : ['supercellmodule.c', 'supercell-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'supercell', '--lower']
)


py.extension_module('supercell',
  [
    'vibes/helpers/supercell/supercell.f90',
    'vibes/helpers/supercell/linalg.f90',
    supercell_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  install_dir: join_paths(py.get_install_dir(), 'vibes/helpers/supercell')
)


# Install all Python files from the vibes directory
install_subdir('vibes',
  install_dir: py.get_install_dir()
)
