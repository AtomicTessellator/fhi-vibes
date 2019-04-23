""" test the harmonic analysis, i.e., mode projection etc. """

import numpy as np
import scipy.linalg as la

from hilde.io import read
from hilde.helpers.numerics import clean_matrix

# from hilde.helpers.supercell import map_indices
from hilde.helpers.lattice_points import (
    map_I_to_iL,
    get_lattice_points,
    get_commensurate_q_points,
)
from hilde.harmonic_analysis import HarmonicAnalysis
from hilde.harmonic_analysis.dynamical_matrix import get_dynamical_matrices
from hilde.harmonic_analysis.displacements import get_U, get_dUdt
from hilde.harmonic_analysis.normal_modes import u_I_to_u_s, u_s_to_u_I, get_A_qst2
from hilde.trajectory import reader
from hilde.tdep.wrapper import parse_tdep_forceconstant

from hilde.konstanten import kB
from hilde.helpers import Timer

from harmonic_md import run


def main():
    timer = Timer()

    primitive = read("geometry.in.primitive")
    supercell = read("geometry.in.supercell")
    force_constants = parse_tdep_forceconstant("infile.forceconstant_remapped")

    masses = supercell.get_masses()

    lattice_points, _ = get_lattice_points(primitive.cell, supercell.cell)
    indeces = map_I_to_iL(primitive, supercell)

    # check if the commensurate q point is correct
    q_points = get_commensurate_q_points(primitive.cell, supercell.cell)
    assert la.norm(q_points[1] - [0.18429553175324292, 0.0, 0.0]) < 1e-14, q_points

    # diagonalize dynamical matrices at commensurate q points

    dyn_matrices = get_dynamical_matrices(
        q_points, primitive, supercell, force_constants
    )

    omegas2, evs = [], []
    for k, dyn_matrix in zip(q_points, dyn_matrices):
        w_2, ev = la.eigh(dyn_matrix)
        omegas2.append(w_2)
        evs.append(ev)

    omegas2 = np.array(omegas2)
    evs = np.array(evs)

    aux_args = {
        "q_points": q_points,
        "lattice_points": lattice_points,
        "eigenvectors": evs,
        "indeces": indeces,
    }

    # check if eigenvectors are orthogonal
    for iq, q in enumerate(q_points):
        e = evs[iq, :, :]
        diff = e.conj().T @ e - np.eye(evs.shape[1])
        assert la.norm(diff) < 1e-14, (q, evs)

    # check if transformation is unitary by probing each mode

    U = np.zeros_like(supercell.positions)
    u_qs = np.zeros((4, 6))

    for (i, j) in np.ndindex(u_qs.shape):
        u_qs *= 0
        u_qs[i, j] = -4
        u_I = u_s_to_u_I(u_qs, **aux_args)

        assert la.norm(u_qs - u_I_to_u_s(u_I, **aux_args)) < 1e-14, (u_qs, u_I)

    # set velocities such that temperature is 100K
    temp = 100
    omegas = omegas2 ** 0.5
    pref = (2 * kB * temp) ** 0.5
    amplitudes = pref / omegas * (omegas.size / (omegas.size - 3)) ** 0.5

    # set acoustic modes to zero
    amplitudes[0, :3] = 0

    const = 1 / (2 * kB) / 4 / 2 / 3
    assert la.norm(const * (amplitudes ** 2 * omegas ** 2).sum() - temp) < 1e-13

    # \dot u = \omega * A
    V = u_s_to_u_I(omegas * amplitudes, **aux_args)

    # mass scaling
    V /= masses[:, None] ** 0.5

    prepared_cell = supercell.copy()
    prepared_cell.set_velocities(V)

    # check that the temperature is as expected
    prepared_temp = prepared_cell.get_temperature()
    assert la.norm(prepared_temp - 2 * temp) / temp < 1e-5, prepared_temp

    # write prepared cell as input for MD and run
    prepared_cell.write("geometry.in", format="aims", velocities=True)

    run(harmonic=True, maxsteps=501, dt=2)

    # read the obtained trajectory and check the average temperature
    traj = reader("trajectory.yaml")

    temperatures = np.array([a.get_temperature() for a in traj])

    assert abs((temperatures.mean() - temp) / temp) < 0.001, temperatures.mean()

    # Test if amplitudes are directly restored

    atoms_displaced = traj[0]

    u_qst = [u_I_to_u_s(get_U(atoms_displaced, supercell), **aux_args)]
    v_qst = [u_I_to_u_s(get_dUdt(atoms_displaced), **aux_args)]

    new_amplitudes = get_A_qst2(u_qst, v_qst, omegas2) ** 0.5
    assert la.norm(new_amplitudes - amplitudes) < 1e-10, (new_amplitudes, amplitudes)

    # check that mode projection preserves kinetic energy
    U_t = [u_I_to_u_s(get_U(atoms, supercell), **aux_args) for atoms in traj]
    V_t = [u_I_to_u_s(get_dUdt(atoms), **aux_args) for atoms in traj]

    const = 1 / kB / 4 / 2 / 3

    for ii in range(100):
        amp = V_t[ii]
        t = const * (amp ** 2).sum()
        assert abs(t - traj[ii].get_temperature()) < 1e-4, t

    # Check that energy in each mode was approx. constant

    a = get_A_qst2(U_t, V_t, omegas2)

    E = 0.5 * omegas[None, :, :] ** 2 * a

    assert abs(E.mean() / kB - temp) / temp < 0.01, E.mean() / kB
    assert E[:, 3:, :].std() / E.mean() < 0.01, E[:, 3:, :].std() / E.mean()

    # compare the high level access via HarmonicAnalysis
    fcs, lps = parse_tdep_forceconstant("infile.forceconstant")

    ha = HarmonicAnalysis(primitive, supercell, fcs, lps)
    A, p, E = ha.project(traj)

    # check that eigenvectors coincide
    for ii, q in enumerate(ha.q_points_frac):
        _, e1 = ha.solve_Dq(q)
        e2 = evs[ii]
        assert la.norm(e1 - e2) < 1e-14, (e1, e2)

    # make sure the mode energies were conserved
    assert abs(E.mean() / kB - temp) / temp < 0.01, E.mean() / kB
    assert E[:, 3:, :].std() / E.mean() < 0.01, E[:, 3:, :].std() / E.mean()

    timer("ran mode projection test successfully")


if __name__ == "__main__":
    main()