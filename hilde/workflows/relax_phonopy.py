''' An example of how to use FireWorks in conjunction with HilDe and remote clusters'''
from argparse import ArgumentParser
import os

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire

from hilde.fireworks_api_adapter.qlaunch_remote import qlaunch_remote
from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.parsers import read_structure
from hilde.phonon_db.phonon_db import connect
from hilde.structure.structure import patoms2dict, calc2dict
from hilde.tasks import fireworks as fw
from hilde.tasks.fireworks import mutate_kgrid
from hilde.tasks.fireworks import energy_diff
from hilde.templates.aims import setup_aims

aims_kgrid_conv_settings = {
    "species_type" : "light",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_iter_limit": 1000
}
aims_relax_settings_light = {
    "species_type" : "light",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "relax_geometry": "trm 1E-2",
    "relax_unit_cell": "full",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}
aims_relax_settings_tight = {
    "species_type" : "tight",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "relax_geometry": "trm 1E-2",
    "relax_unit_cell": "full",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}
aims_force_settings = {
    "species_type" : "tight",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}

def gen_kgrid_conv_fw(atoms,
                      workdir,
                      atoms_spec,
                      final_atoms_spec,
                      name="k_grid_conv",
                      spec_qad={},
                      calc=None,
                      calc_settings=aims_kgrid_conv_settings):
    if calc is None:
        calc = setup_aims(calc_settings)
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
                            "kwargs": {"criteria": 1e-3}}))
    return Firework(task_list, name=name)

def gen_relax_fw(atoms,
                 db_name,
                 workdir,
                 in_atoms_spec,
                 out_atoms_spec,
                 calc=None,
                 up_calc_from_db=None,
                 name="relax",
                 spec_qad={},
                 from_db=False,
                 calc_settings=aims_relax_settings_light):
    if calc is None:
        calc = setup_aims(calc_settings)
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
    if up_calc_from_db:
        inputs_calc.append("calculator")
        task_list.append(PyTask({"func": fw.mod_calc.name,
                                 "args": [up_calc_from_db[0], calc2dict(calc)],
                                 "inputs": [up_calc_from_db[0]],
                                 "kwargs": {"spec_key": up_calc_from_db[0]}}))
        for key in up_calc_from_db[1:]:
            task_list.append(PyTask({"func": fw.mod_calc.name,
                                     "args": [key],
                                     "inputs": ["calculator",key],
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
                             "kwargs": {"fw_name": name, "original_atoms_hash": atoms_hash}}))
    return Firework(task_list, name=name)

def gen_initialize_phonopy_fw(atoms,
                              smatrix,
                              workdir,
                              atoms_spec=None,
                              calc=None,
                              symprec=1e-5,
                              up_calc_from_db=None,
                              name="init_phono",
                              spec_qad={},
                              from_db=False,
                              calc_settings=aims_force_settings):
    if calc is None:
        calc = setup_aims(calc_settings)
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
                                     "inputs": ["calculator",key]}))
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
                           workdir,
                           atoms_spec=None,
                           symprec=1e-5,
                           name="analyze_phono",
                           from_db=False):
    atoms_hash, calc_hash = hash_atoms(atoms)
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

    task_list = []
    task_list.append(PyTask({"func": fw.calc_phonopy_force_constants.name,
                             "args": args_fc,
                             "inputs": inputs_fc}))
    task_list.append(PyTask({"func": fw.add_phonon_to_db.name,
                             "args": args_to_db,
                             "inputs": inputs_to_db,
                             "kwargs": {"symprec": symprec, "fw_name": name, "original_atoms_hash": atoms_hash}}))
    return Firework(task_list, name=name)

def gen_relax_phonopy_wf(geo_in_file,
                         db_name_remote,
                         db_name_local,
                         name,
                         workdir,
                         atoms_spec,
                         smatrix,
                         symprec=1e-5,
                         kgrid_conv=aims_kgrid_conv_settings,
                         relax_light=aims_relax_settings_light,
                         relax_tight=aims_relax_settings_tight,
                         force_calc=aims_force_settings,
                         spec_qad_kgrid={},
                         spec_qad_relax={},
                         spec_qad_forces={}):
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
                                 workdir + "/force_calcs/",
                                 "tight_relax_atoms",
                                 symprec=symprec,
                                 name=f"analyze_phono_{name}",
                                 from_db=True)
    workflow = Workflow([fw1, fw2, fw3, fw4,fw5], {fw1:[fw2], fw2:[fw3], fw3:[fw4], fw4:[fw5]},
                        name=name)
    return workflow
