'''Generates a Workflow to relax a structure and calculate its harmonic force constants'''
from fireworks import Firework, PyTask, Workflow

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import setup_workdir
from hilde.parsers import read_structure
from hilde.structure.structure import patoms2dict, calc2dict
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
                      criteria=1e-3):
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
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings)
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
                 calc_settings=None):
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
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings)
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
                             "inputs": inputs_calc,
                             "spec": spec_qad}))
    task_list.append(PyTask({"func": fw.get_relaxed_structure.name,
                             "args": [str(workdir/"geometry.in.next_step"), out_atoms_spec],
                             "inputs": [in_atoms_spec]}))
    task_list.append(PyTask({"func": fw.add_phonon_to_db.name,
                             "args": [db_name],
                             "inputs": [in_atoms_spec, out_atoms_spec],
                             "kwargs": kwargs_to_db}))
    return Firework(task_list, name=name)

def gen_initialize_phonopy_fw(atoms,
                              smatrix,
                              workdir,
                              atoms_spec=None,
                              calc=None,
                              symprec=1e-5,
                              up_calc_from_db=None,
                              name="init_phono",
                              spec_qad=None,
                              from_db=False,
                              calc_settings=None):
    '''
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
    '''
    if calc is None:
        calc = setup_aims(custom_settings=calc_settings)
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
        task_list.append(PyTask({"func": fw.mod_calc.name,
                                 "args": [up_calc_from_db[0], calc2dict(calc)],
                                 "inputs": [up_calc_from_db[0]]}))
        for key in up_calc_from_db[1:]:
            task_list.append(PyTask({"func": fw.mod_calc.name,
                                     "args": [key],
                                     "inputs": ["calculator", key]}))
    kwargs_calc = {"spec_qad": spec_qad, "out_spec": "phonon_calcs"}
    if atoms_spec:
        task_list.append(PyTask({"func": fw.transfer_spec.name,
                                 "args": [atoms_spec],
                                 "inputs": [atoms_spec]}))
    task_list.append(PyTask({"func": fw.initialize_phonopy.name,
                             "args": args_phono,
                             "inputs": inputs_phono,
                             "kwargs": {"symprec": symprec}}))
    task_list.append(PyTask({"func": fw.calculate_multiple.name,
                             "inputs": inputs_calc,
                             "kwargs": kwargs_calc}))
    return Firework(task_list, name=name)

def gen_analyze_phonopy_fw(atoms,
                           db_name,
                           smatrix,
                           atoms_spec=None,
                           db_label="phonons",
                           symprec=1e-5,
                           name="analyze_phono",
                           from_db=False):
    '''
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
    '''
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
    kwargs_to_db = {"symprec": symprec,
                    "fw_name": name,
                    "original_atoms_hash": atoms_hash,
                    "calc_type": db_label}
    task_list = []
    task_list.append(PyTask({"func": fw.calc_phonopy_force_constants.name,
                             "args": args_fc,
                             "inputs": inputs_fc}))
    task_list.append(PyTask({"func": fw.add_phonon_to_db.name,
                             "args": args_to_db,
                             "inputs": inputs_to_db,
                             "kwargs": kwargs_to_db}))
    return Firework(task_list, name=name)

def gen_relax_phonopy_wf(geo_in_file,
                         db_name_remote,
                         db_name_local,
                         name,
                         workdir,
                         smatrix,
                         symprec=1e-5,
                         kgrid_conv=None,
                         relax_light=None,
                         relax_tight=None,
                         force_calc=None,
                         spec_qad_kgrid=None,
                         spec_qad_relax=None,
                         spec_qad_forces=None):
    '''
    Creates a Workflow to relax a structure and find its phonon properties
    Args:
        geo_in_files: str
            geometry.in file
        db_name_remote: str
            The database path on the remote machine
        db_name_local: str
            The database path on the local machine
        name: str
            The name of the Workflow
        workdir: str
            The base work directory to perform the electronic structure calculations
        smatrix: np.ndarray(int64)
            The supercell matrix for the phonon calculations
        symprec: float
            precision used to determine the symmetry in phonopy
        kgrid_conv: dict
            A dictionary of settings used for the kgrid convergence calculations
        relax_light: dict
            A dictionary of settings used for the light relaxation calculations
        relax_tight: dict
            A dictionary of settings used for the tight relaxation calculations
        force_calc: dict
            A dictionary of settings used for the force calculations
        spec_qad_kgrid: dict
            A dictionary describing queue parameters for the job (only changes to the qadapter
            file) for the kgrid convergence calculations
        spec_qad_relax: dict
            A dictionary describing queue parameters for the job (only changes to the qadapter
            file) for the relaxation calculations
        spec_qad_forces: dict
            A dictionary describing queue parameters for the job (only changes to the qadapter
            file) for the force calculations
    Returns: Workflow
        A Workflow that will relax a structure and calculated its phonons
    '''
    if not spec_qad_kgrid:
        spec_qad_kgrid = {}
    if not spec_qad_relax:
        spec_qad_relax = {}
    if not spec_qad_forces:
        spec_qad_forces = {}
    atoms = read_structure(geo_in_file)
    fw1 = gen_kgrid_conv_fw(atoms,
                            workdir + "/kgrid_conv",
                            "atoms_relax",
                            "kgrid_atoms",
                            name=f"k_grid_conv_{name}",
                            spec_qad=spec_qad_kgrid,
                            calc_settings=kgrid_conv)
    fw2 = gen_relax_fw(atoms,
                       db_name_remote,
                       workdir + "/light_relax/",
                       "kgrid_atoms",
                       "light_relax_atoms",
                       db_label="light_relax",
                       up_calc_from_db=["k_grid"],
                       name=f"light_relax_{name}",
                       spec_qad=spec_qad_relax,
                       from_db=False,
                       calc_settings=relax_light)
    fw3 = gen_relax_fw(atoms,
                       db_name_remote,
                       workdir + "/tight_relax/",
                       "light_relax_atoms",
                       "tight_relax_atoms",
                       db_label="tight_relax",
                       up_calc_from_db=["k_grid"],
                       name=f"tight_relax_{name}",
                       spec_qad=spec_qad_relax,
                       from_db=True,
                       calc_settings=relax_tight)
    fw4 = gen_initialize_phonopy_fw(atoms,
                                    smatrix,
                                    workdir + "/force_calcs/",
                                    "tight_relax_atoms",
                                    symprec=symprec,
                                    up_calc_from_db=["k_grid"],
                                    name=f"init_phono_{name}",
                                    spec_qad=spec_qad_forces,
                                    from_db=True,
                                    calc_settings=force_calc)
    fw5 = gen_analyze_phonopy_fw(atoms,
                                 db_name_local,
                                 smatrix,
                                 "tight_relax_atoms",
                                 db_label="phonons",
                                 symprec=symprec,
                                 name=f"analyze_phono_{name}",
                                 from_db=True)
    workflow = Workflow([fw1, fw2, fw3, fw4, fw5], {fw1:[fw2], fw2:[fw3], fw3:[fw4], fw4:[fw5]},
                        name=name)
    return workflow
