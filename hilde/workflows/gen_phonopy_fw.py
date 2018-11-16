"""Generates a Workflow to relax a structure and calculate its harmonic force constants"""
from fireworks import Firework, PyTask, Workflow

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import setup_workdir
from hilde.parsers import read_structure
from hilde.helpers.input_exchange import patoms2dict, calc2dict
from hilde.tasks import fireworks as fw
from hilde.tasks.fireworks import mutate_kgrid
from hilde.templates.aims import setup_aims


def gen_initialize_phonopy_fw(
    atoms,
    smatrix,
    workdir,
    atoms_spec=None,
    calc=None,
    symprec=1e-5,
    up_calc_from_db=None,
    name="init_phono",
    spec_qad=None,
    from_db=False,
    calc_settings=None,
    hilde_cfg="../../hilde.cfg",
):
    """
    Generates a FireWork to initialize the phonon calculations
    Args:
        atoms: pAtoms
            Representation of the initial structure to be relaxed
        smatrix: np.ndarray(int64)
            The supercell matrix for the phonon calculations
        workdir: str
            Path where the calculations should be performed in
        atoms_spec: str
            spec for the FireWorks database (launchpad) where the input atoms should be obtained
        calc: ase calculator
            The calculator that will be modified for the convergence
        symprec: float
            precision used to determine the symmetry in phonopy
        up_calc_from_db: list of strs
            A list of spec keys to be used to modify the calculator
        name: str:
            Name for the FireWork
        spec_qad: dict
            A dictionary describing queue parameters for the job (only changes to the qadapter
            file should be included)
        from_db: bool
            if True use in_atoms_spec to get the atoms structure to be relaxed
        calc_settings: dictionary
            A dictionary of custom settings that should be used to make an aims calculator
    Returns: FireWork
        A FireWork that will initialize a phonopy calculation and add the force calculations to
        the Workflow
    """
    if calc_settings is None:
        calc_settings = {}
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings, config_file=hilde_cfg)
    if not spec_qad:
        spec_qad = {}
    atoms.set_calculator(calc)
    workdir = setup_workdir(atoms, workdir, False)
    atoms = patoms2dict(atoms)

    args_phono = [smatrix, workdir]
    inputs_phono = []
    inputs_calc = ["workdirs", "atom_dicts"]
    task_list = []
    if from_db:
        inputs_phono.append(atoms_spec)
    else:
        args_phono.append(atoms)
    if up_calc_from_db:
        inputs_calc.append("calculator")
        task_list.append(
            PyTask(
                {
                    "func": fw.mod_calc.name,
                    "args": [up_calc_from_db[0], calc2dict(calc)],
                    "inputs": [up_calc_from_db[0]],
                }
            )
        )
        for key in up_calc_from_db[1:]:
            task_list.append(
                PyTask({"func": fw.mod_calc.name, "args": [key], "inputs": ["calculator", key]})
            )
    kwargs_calc = {"spec_qad": spec_qad, "out_spec": "phonon_calcs"}
    if atoms_spec:
        task_list.append(
            PyTask({"func": fw.transfer_spec.name, "args": [atoms_spec], "inputs": [atoms_spec]})
        )

    task_list.append(
        PyTask(
            {
                "func": fw.initialize_phonopy.name,
                "args": args_phono,
                "inputs": inputs_phono,
                "kwargs": {"symprec": symprec},
            }
        )
    )
    task_list.append(
        PyTask({"func": fw.calculate_multiple.name, "inputs": inputs_calc, "kwargs": kwargs_calc})
    )
    return Firework(task_list, name=name)


def gen_analyze_phonopy_fw(
    atoms,
    db_name,
    smatrix,
    atoms_spec=None,
    db_label="phonons",
    symprec=1e-5,
    name="analyze_phono",
    from_db=False,
):
    """
    Generates a FireWork to calculate the second order force constants
    Args:
        atoms: pAtoms
            Representation of the initial structure to be relaxed
        db_name: str
            Path to the database where to store the results
        smatrix: np.ndarray(int64)
            The supercell matrix for the phonon calculations
        atoms_spec: str
            spec for the FireWorks database (launchpad) where the input atoms should be obtained
        db_label: str
            Label describing the calculation type
        symprec: float
            precision used to determine the symmetry in phonopy
        name: str:
            Name for the FireWork
        from_db: bool
            If True get the initial structure from the database
    Returns: FireWork
        A FireWork that will initialize a phonopy calculation and add the force calculations to
        the Workflow
    """
    atoms_hash, _ = hash_atoms(atoms)
    atoms = patoms2dict(atoms)
    args_fc = [smatrix]
    args_to_db = [db_name]
    if from_db:
        inputs_fc = [atoms_spec, "phonon_calcs"]
        inputs_to_db = [atoms_spec, "phonon_dict"]
    else:
        args_fc.append(atoms)
        args_to_db.append(atoms)
        inputs_fc = ["phonon_calcs"]
        inputs_to_db = ["phonon_dict"]

    kwargs_to_db = {
        "symprec": symprec,
        "fw_name": name,
        "original_atoms_hash": atoms_hash,
        "calc_type": db_label,
    }
    task_list = []
    task_list.append(
        PyTask({"func": fw.calc_phonopy_force_constants.name, "args": args_fc, "inputs": inputs_fc})
    )
    task_list.append(
        PyTask(
            {
                "func": fw.add_phonon_to_db.name,
                "args": args_to_db,
                "inputs": inputs_to_db,
                "kwargs": kwargs_to_db,
            }
        )
    )
    return Firework(task_list, name=name)
