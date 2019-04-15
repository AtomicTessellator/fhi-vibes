import numpy as np

from hilde.helpers.pickle import pread
from hilde.phonopy.postprocess import postprocess
from hilde.phonopy.wrapper import get_dos
from hilde.phonopy.wrapper import get_bandstructure

ph_file = "phonopy/phonon.pick.gz"
ph = pread(ph_file)
total_dos = get_dos(ph, write=True)
dos = np.loadtxt("total_dos.dat")
assert np.abs(np.max(dos[:,0] - total_dos["frequency_points"])) < 1e-13
assert np.abs(np.max(dos[:,1] - total_dos["total_dos"])) < 1e-10

phonon_bs = get_bandstructure(ph)
