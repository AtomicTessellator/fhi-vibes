# flake8: noqa
from ase import units
from numpy import pi


v_unit = units.Ang / (1000.0 * units.fs)

# Symmetry
symprec = 1e-5

# io
n_db_digits = 14
n_geom_digits = 12
n_yaml_digits = 14

# maths
perfect_fill = 0.523598775598299
vol_sphere = 4.0 * pi / 3.0

# physical constants
AMU = units._amu  # = 1.66053904e-27  # [kg]
LIGHT = units._c  # = 299792458  # [m / s]
PLANCK_CONSTANT = units._hplanck  # = 6.62607015e-34  # [J s]
BOLTZMANN = units._k  # = 1.38064852e-23  # [J / K]
AVOGADRO = units._Nav  # = 6.02214076e23  # [1]
ALPHA = units.alpha  # = 1 / 137.035999046  # [1]

# Mathematical constants
PI = pi

# Conversion factors
AA = 1e-10  # [m]
PICO = 1e-12  # [s]
FEMTO = 1e-3 * PICO  # [s]
EV = units._e  # = 1.60217733e-19  # [J]
THZ = 1 / PICO  # [1/s]

# Atomic units
ELECTRON_MASS = units._me  # = 5.48579909070e-4 * AMU  # [kg]
HARTREE = units.Hartree * EV  # = 27.21138602 * EV  # [J]
BOHR = units.Bohr * AA  # = 0.52917721092 * AA  # [m]
HBAR = units._hbar  # = PLANCK_CONSTANT / 2 / PI  # [J s]

# phonons related
omega_to_THz = (EV / AA ** 2 / AMU) ** 0.5 / THZ / 2 / PI  # 15.633302 THz
THz_to_cm = THZ / LIGHT / 100  # 33.3564 [1/cm]
omega_to_cm = omega_to_THz * THz_to_cm
gv_to_AA_fs = omega_to_THz / 1000  # group velocity to AA/fs

amu_AA_THz_to_eV = AMU * AA ** 2 * THZ ** 2 / EV

EvTokJmol = EV / 1000 * AVOGADRO  # [kJ/mol] 96.4853910
kJmolToEv = 1 / EvTokJmol
THzToEv = PLANCK_CONSTANT * 1e12 / EV  # [eV]
