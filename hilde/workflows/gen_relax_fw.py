'''Generates a Workflow to relax a structure and calculate its harmonic force constants'''
from fireworks import Firework, PyTask, Workflow

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import setup_workdir
from hilde.parsers import read_structure
from hilde.helpers.input_exchange import patoms2dict, calc2dict
from hilde.tasks import fireworks as fw
from hilde.tasks.fireworks import mutate_kgrid
from hilde.templates.aims import setup_aims

def gen_kgrid_conv_fw(atoms,
                      workdir,
                      atoms_spec,
                      final_atoms_spec,
                      name="k_grid_conv",
                      spec_qad=None,
                      calc=None,
                      calc_settings=None,
                      criteria=1e-3,
                      hilde_cfg="../../hilde.cfg"):
    '''
    Generates a k_grid convergence FireWork
    Args:
        atoms: pAtoms
            Representation of the structure we want to converge with
        workdir: str
            Path where the calculations should be performed in
        atoms_spec: str
            spec for the FireWorks database (launchpad) where the calculated atoms should be stored
        final_atoms_spec: str
            spec for the FireWorks database where the converged atoms will be stored
        name: str:
            Name for the FireWork
        spec_qad: dict
            A dictionary describing queue parameters for the job (only changes to the qadapter
            file should be included)
        calc: ase calculator
            The calculator that will be modified for the convergence
        calc_settings: dictionary
            A dictionary of custom settings that should be used to make an aims calculator
        criteria: float
            Convergence threshold
    Returns: FireWork
        A FireWork that will do a kgrid convergence for your system
    '''
    if calc_settings is None:
        calc_settings = {}
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings, config_file=hilde_cfg)
    if not spec_qad:
        spec_qad = {}
    atoms.set_calculator(calc)
    workdir = setup_workdir(atoms, workdir, False)
    workdirs = [str(workdir/f'{0:05d}'), str(workdir/f'{1:05d}')]
    atoms_dicts = [patoms2dict(atoms), mutate_kgrid(patoms2dict(atoms))[0]]
    task_list = []
    task_list.append(PyTask({"func": fw.calculate.name,
                             "args": [workdirs[0], atoms_spec, atoms_dicts[0]],
                             "spec": spec_qad}))
    task_list.append(PyTask({"func": fw.calculate.name,
                             "args": [workdirs[1], atoms_spec, atoms_dicts[1]],
                             "spec": spec_qad}))
    conv_args = ["k_grid",
                 atoms_spec,
                 final_atoms_spec,
                 fw.energy_diff.name,
                 fw.mutate_kgrid.name,
                 workdirs[1]]
    task_list.append(PyTask({"func": fw.check_convergence.name,
                             "args": conv_args,
                             "inputs": [atoms_spec],
                             "kwargs": {"criteria": criteria}}))
    return Firework(task_list, name=name)

def gen_relax_fw(atoms,
                 db_name,
                 workdir,
                 in_atoms_spec,
                 out_atoms_spec,
                 db_label='relax',
                 calc=None,
                 up_calc_from_db=None,
                 name="relax",
                 spec_qad=None,
                 from_db=False,
                 calc_settings=None,
                 hilde_cfg="../../hilde.cfg"):
    '''
    Generates a relaxation FireWork
    Args:
        atoms: pAtoms
            Representation of the initial structure to be relaxed
        db_name: str
            Path to the database where to store the atoms
        workdir: str
            Path where the calculations should be performed in
        in_atoms_spec: str
            spec for the FireWorks database (launchpad) where the input atoms should be obtained
        out_atoms_spec: str
            spec for the FireWorks database where the relaxed atoms will be stored
        db_label: str
            Label describing the calculation type
        calc: ase calculator
            The calculator that will be modified for the convergence
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
        A FireWork that will relax a structure with a given set of parameters
    '''
    if calc_settings is None:
        calc_settings = {}
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings, config_file=hilde_cfg)
    if not spec_qad:
        spec_qad = {}
    atoms.set_calculator(calc)
    workdir = setup_workdir(atoms, workdir, False)
    atoms_hash, calc_hash = hash_atoms(atoms)
    atoms = patoms2dict(atoms)
    task_list = []
    args_calc = [workdir, "temporary_atoms_push"]
    inputs_calc = []
    if from_db:
        inputs_calc.append(in_atoms_spec)
    else:
        args_calc.append(atoms)
    kwargs_to_db = {"fw_name": name, "original_atoms_hash": atoms_hash,"calc_type": db_label}
    if up_calc_from_db:
        inputs_calc.append("calculator")
        task_list.append(PyTask({"func": fw.mod_calc.name,
                                 "args": [up_calc_from_db[0], calc2dict(calc)],
                                 "inputs": [up_calc_from_db[0]],
                                 "kwargs": {"spec_key": up_calc_from_db[0]}}))
        for key in up_calc_from_db[1:]:
            task_list.append(PyTask({"func": fw.mod_calc.name,
                                     "args": [key],
                                     "inputs": ["calculator", key],
                                     "kwargs": {"spec_key": key}}))

    task_list.append(PyTask({"func": fw.calculate.name,
                             "args": args_calc,
                             "inputs": inputs_calc}))
    task_list.append(PyTask({"func": fw.get_relaxed_structure.name,
                             "args": [str(workdir/"geometry.in.next_step"), out_atoms_spec],
                             "inputs": [in_atoms_spec]}))
    task_list.append(PyTask({"func": fw.add_phonon_to_db.name,
                             "args": [db_name],
                             "inputs": [in_atoms_spec, out_atoms_spec],
                             "kwargs": kwargs_to_db}))
    return Firework(task_list, name=name, spec=spec_qad)