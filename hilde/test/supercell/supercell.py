from hilde.structure import read_aims
from hilde.helpers.geometry import get_cubicness
from hilde.helpers.linalg import clean_matrix
import numpy as np
gao = read_aims('gan.in')

gao.inform()
print()

for nn in [10]: #range(10, 500, 10):
    smatrix = gao.find_cubic_cell(Ntarget=nn, deviation=.2, verbose=0)
    scell = clean_matrix(smatrix @ gao.cell)
    cness = get_cubicness(scell)
    Nreached = np.linalg.det(smatrix)*gao.n_atoms
    print(f'N target, N reached, cubicness: {nn:4d} {Nreached:7.2f} {cness:5.3f} ({cness**3:5.3f})')
    print(f'Smatrix: {smatrix.flatten()}')