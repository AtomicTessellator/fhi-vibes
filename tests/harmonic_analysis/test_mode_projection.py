""" test the harmonic analysis, i.e., mode projection etc. """
from pathlib import Path

import numpy as np
import scipy.linalg as la
from ase import units
from ase.md.verlet import VelocityVerlet

from vibes.ase.calculators.fc import FCCalculator
from vibes.harmonic_analysis import HarmonicAnalysis
from vibes.harmonic_analysis.dynamical_matrix import get_dynamical_matrices
from vibes.harmonic_analysis.normal_modes import get_A_qst2, projector, u_s_to_u_I
from vibes.helpers import Timer, progressbar
from vibes.helpers.displacements import get_dUdt, get_U
from vibes.helpers.lattice_points import get_lattice_points, map_I_to_iL
from vibes.helpers.supercell import get_commensurate_q_points
from vibes.io import read
from vibes.konstanten import kB
from vibes.molecular_dynamics.utils import MDLogger
from vibes.tdep.wrapper import (
    parse_tdep_forceconstant,
    parse_tdep_remapped_forceconstant,
)
from vibes.trajectory import reader

parent = Path(__file__).parent


def test_all():
    timer = Timer()

    primitive = read(parent / "geometry.in.primitive")
    supercell = read(parent / "geometry.in.supercell")
    force_constants = parse_tdep_forceconstant(
        fc_file=parent / "infile.forceconstant",
        primitive=parent / "geometry.in.primitive",
        supercell=parent / "geometry.in.supercell",
        two_dim=True,
        format="aims",
    )

    masses = supercell.get_masses()

    lattice_points, _ = get_lattice_points(primitive.cell, supercell.cell)
    indeces, _ = map_I_to_iL(primitive, supercell)

    # check if the commensurate q point is correct
    q_points = get_commensurate_q_points(primitive.cell, supercell.cell)
    # assert la.norm(q_points[1] - [0.184_295_531_753_242_92, 0.0, 0.0]) < 1e-14, q_points

    # diagonalize dynamical matrices at commensurate q points

    dyn_matrices = get_dynamical_matrices(
        q_points, primitive, supercell, force_constants
    )

    omegas2, evs = [], []
    for _, dyn_matrix in zip(q_points, dyn_matrices):
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

    P = projector(**aux_args)

    # check if eigenvectors are orthogonal
    for iq, q in enumerate(q_points):
        e = evs[iq, :, :]
        diff = e.conj().T @ e - np.eye(evs.shape[1])
        assert la.norm(diff) < 1e-14, (q, evs)

    # check if transformation is unitary by probing each mode
    u_qs = np.zeros((4, 6))

    for (i, j) in np.ndindex(u_qs.shape):
        u_qs *= 0
        u_qs[i, j] = -4
        u_I = u_s_to_u_I(u_qs, **aux_args)

        diff = la.norm(u_qs - P @ u_I.flatten())
        assert diff < 1e-14, (u_qs, u_I)

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

    _run_md(harmonic=True, maxsteps=501, dt=2, sample="geometry.in")

    # read the obtained trajectory and check the average temperature
    traj = reader(parent / "trajectory.son", verbose=False)

    temperatures = np.array([a.get_temperature() for a in traj])

    assert abs((temperatures.mean() - temp) / temp) < 0.001, temperatures.mean()

    # Test if amplitudes are directly restored

    atoms_displaced = traj[0]

    u_qst = [P @ get_U(atoms_displaced, supercell).flatten()]
    v_qst = [P @ get_dUdt(atoms_displaced).flatten()]

    new_amplitudes = get_A_qst2(u_qst, v_qst, omegas2) ** 0.5
    assert la.norm(new_amplitudes - amplitudes) < 1e-10, (new_amplitudes, amplitudes)

    # check that mode projection preserves kinetic energy
    U_t = [P @ get_U(atoms, supercell).flatten() for atoms in traj]
    V_t = [P @ get_dUdt(atoms).flatten() for atoms in traj]

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
    fcs, lps = parse_tdep_remapped_forceconstant(parent / "infile.forceconstant")

    ha = HarmonicAnalysis(primitive, supercell, fcs, lps)
    _, _, E = ha.project(traj)

    # check that eigenvectors coincide
    for ii, q in enumerate(ha.q_points_frac):
        _, e1 = ha.solve_Dq(q)
        e2 = evs[ii]
        assert la.norm(e1 - e2) < 1e-14, (e1, e2)

    # make sure the mode energies were conserved
    assert abs(E.mean() / kB - temp) / temp < 0.01, E.mean() / kB
    assert E[:, 3:, :].std() / E.mean() < 0.01, E[:, 3:, :].std() / E.mean()

    timer("ran mode projection test successfully")


def _run_md(
    maxsteps=1001,
    dt=1,
    harmonic=True,
    sample=parent / "geometry.in.supercell.300K",
    primitive=parent / "geometry.in.primitive",
    supercell=parent / "geometry.in.supercell",
    fc_file=parent / "infile.forceconstant",
    trajectory_file=parent / "trajectory.son",
):
    """ run Verlet MD, harmonic or force field """
    trajectory_file = trajectory_file
    atoms = read(sample)

    force_constants = parse_tdep_forceconstant(
        fc_file=fc_file,
        primitive=primitive,
        supercell=supercell,
        two_dim=True,
        format="aims",
    )
    # force_constants.resize(2 * (3 * len(supercell),))

    supercell = read(supercell)
    if harmonic is True:
        calculator = FCCalculator(supercell, force_constants)
    else:
        raise RuntimeError("FIXME")

    # generic md settings
    settings = {"atoms": atoms, "timestep": dt * units.fs}
    metadata = {"MD": {"fs": units.fs, "dt": dt}}

    md = VelocityVerlet(**settings)

    logger = MDLogger(atoms, trajectory_file, metadata=metadata, overwrite=True)

    atoms.calc = calculator
    for _ in progressbar(range(maxsteps)):
        logger(atoms, info={"nsteps": md.nsteps, "dt": md.dt})
        md.run(1)


if __name__ == "__main__":
    import ase

    if ase.__version__ >= "3.18":
        test_all()
    else:
        print(f"ase version: {ase.__version__}, skip for now")
