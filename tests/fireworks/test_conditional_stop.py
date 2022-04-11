import numpy as np
import sys

from ase.calculators.emt import EMT
from ase.io.aims import read_aims
from fireworks import FWAction
from pathlib import Path

from vibes.fireworks.tasks.fw_out.check_conditionals import (
    check_stop_conditional,
    check_condition,
    check_ojbect,
    run_all_checks,
)
from vibes.phonopy.postprocess import postprocess
from vibes.trajectory.io import reader

parent = Path(__file__).parent


def test_conditions():
    assert check_condition(2.0, ["", "gt", 1.0])
    assert check_condition(2.0, ["", "ge", 1.0])
    assert check_condition(1.0, ["", "ge", 1.0])
    assert check_condition(1.0, ["", "le", 1.0])
    assert check_condition(0.0, ["", "le", 1.0])
    assert check_condition(0.0, ["", "lt", 1.0])
    assert check_condition(0.0, ["", "ne", 1.0])
    assert check_condition(1.0, ["", "eq", 1.0])

    assert check_condition([0, 1, 2], ["max", "eq", 2.0])
    assert check_condition([0, 1, 2], ["min", "eq", 0.0])
    assert check_condition([0, 1, 2], ["mean", "eq", 1.0])
    assert check_condition([0, 1, 2], ["median", "eq", 1.0])
    assert check_condition([2, 2, 2], ["std", "eq", 0.0])


def test_user_defined_conditions():
    atoms = read_aims(f"{parent}/geometry.in")
    sys.path.append(f"{parent}/conditional_functions/")

    assert check_stop_conditional("cond_test_va.cond_test_atomic_volume", atoms)
    atoms.set_cell(100 * np.eye(3))
    assert not check_stop_conditional("cond_test_va.cond_test_atomic_volume", atoms)


def test_aims_stop():
    atoms = read_aims(f"{parent}/geometry.in")
    atoms.set_calculator(EMT())
    condition_list = [
        ["volume", "", "lt", 0.0],
        ["energy", "", "lt", 0.0],
        ["free_energy", "", "lt", 0.0],
        ["lattice_parameters", "max", "lt", 6.79],
        ["lattice_angles", "max", "le", 100],
    ]
    assert not check_ojbect(atoms, condition_list)

    condition_list[0][2] = "gt"
    assert check_ojbect(atoms, condition_list)

    condition_list[0][2] = "lt"
    condition_list[1][2] = "gt"
    assert check_ojbect(atoms, condition_list)

    condition_list[1][2] = "lt"
    condition_list[2][2] = "gt"
    assert check_ojbect(atoms, condition_list)

    condition_list[2][2] = "lt"
    condition_list[3][2] = "gt"
    assert check_ojbect(atoms, condition_list)

    condition_list[3][2] = "lt"
    condition_list[4][2] = "gt"
    assert check_ojbect(atoms, condition_list)


def test_phonopy_stop():
    phonon = postprocess(f"{parent}/conditional_check_trajs/phonon.son")
    condition_list = [
        ["cv", "", "gt", 40.0, 300],
        ["free_energy", "", "lt", 0.0, 300],
        ["entropy", "", "gt", 45.0, 300],
        ["frequencies", "min", "gt", 0.001, [0, 0, 0]],
        ["theta_D", "", "gt", 400],
        ["theta_D_infty", "", "gt", 700],
        ["theta_P", "", "gt", 500],
    ]
    assert not check_ojbect(phonon, condition_list)

    condition_list[0][4] = 600
    assert check_ojbect(phonon, condition_list)

    condition_list[0][4] = 300
    condition_list[1][4] = 600
    assert check_ojbect(phonon, condition_list)

    condition_list[1][4] = 300
    condition_list[2][4] = 600
    assert check_ojbect(phonon, condition_list)

    condition_list[2][4] = 300
    condition_list[3][4] = [0.5, 0.0, 0.0]
    assert check_ojbect(phonon, condition_list)

    condition_list[3][4] = [0.0, 0.0, 0.0]
    condition_list[4][2] = "lt"
    assert check_ojbect(phonon, condition_list)

    condition_list[4][2] = "gt"
    condition_list[5][2] = "lt"
    assert check_ojbect(phonon, condition_list)

    condition_list[5][2] = "gt"
    condition_list[6][2] = "lt"
    assert check_ojbect(phonon, condition_list)


def test_stat_samp_stop():
    trajectory = reader(f"{parent}/conditional_check_trajs/stat_sampling.son")
    condition_list = [["sigma", "", "gt", 0.13]]
    assert not check_ojbect(trajectory, condition_list)

    condition_list[0][2] = "lt"
    assert check_ojbect(trajectory, condition_list)


def test_all_checks():
    atoms = read_aims(f"{parent}/geometry.in")
    atoms.set_calculator(EMT())
    stop_if = {
        "condition_list": [
            ["volume", "", "lt", 0.0],
            ["energy", "", "lt", 0.0],
            ["free_energy", "", "lt", 0.0],
            ["lattice_parameters", "max", "lt", 6.79],
            ["lattice_angles", "max", "le", 100],
        ],
        "external_functions": ["cond_test_va.cond_test_atomic_volume"],
    }
    assert isinstance(run_all_checks(atoms, stop_if), FWAction)

    stop_if["external_functions"] = []
    assert run_all_checks(atoms, stop_if) is None

    del stop_if["external_functions"]
    assert run_all_checks(atoms, stop_if) is None

    stop_if["external_functions"] = []
    stop_if["condition_list"][0][2] = "gt"
    assert isinstance(run_all_checks(atoms, stop_if), FWAction)

    stop_if["condition_list"] = []
    assert run_all_checks(atoms, stop_if) is None

    del stop_if["condition_list"]
    assert run_all_checks(atoms, stop_if) is None

    stop_if["external_functions"] = ["cond_test_va.cond_test_atomic_volume"]
    assert isinstance(run_all_checks(atoms, stop_if), FWAction)


if __name__ == "__main__":
    test_conditions()
    test_user_defined_conditions()
    test_aims_stop()
    test_phonopy_stop()
    test_stat_samp_stop()
    test_all_checks()
