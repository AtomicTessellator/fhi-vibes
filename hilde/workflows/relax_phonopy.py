'''Generates a Workflow to relax a structure and calculate its harmonic force constants'''
from fireworks import Firework, PyTask, Workflow

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import setup_workdir
from hilde.parsers import read_structure
from hilde.helpers.input_exchange import patoms2dict, calc2dict
from hilde.tasks import fireworks as fw
from hilde.tasks.fireworks import mutate_kgrid
from hilde.templates.aims import setup_aims
from hilde.workflows.gen_relax_fw import gen_kgrid_conv_fw, gen_relax_fw
from hilde.workflows.gen_phonopy_fw import gen_initialize_phonopy_fw, gen_analyze_phonopy_fw

def gen_relax_phonopy_wf(
    geo_in_file,
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
    spec_qad_forces=None,
    hilde_cfg="../../hilde.cfg",
):
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
    fw1 = gen_kgrid_conv_fw(
        atoms,
        workdir + "/kgrid_conv",
        "atoms_relax",
        "kgrid_atoms",
        name=f"k_grid_conv_{name}",
        spec_qad=spec_qad_kgrid,
        calc_settings=kgrid_conv,
        hilde_cfg=hilde_cfg
    )
    fw2 = gen_relax_fw(
        atoms,
        db_name_remote,
        workdir + "/light_relax/",
        "kgrid_atoms",
        "light_relax_atoms",
        db_label="light_relax",
        up_calc_from_db=["k_grid"],
        name=f"light_relax_{name}",
        spec_qad=spec_qad_relax,
        from_db=False,
        calc_settings=relax_light,
        hilde_cfg=hilde_cfg
    )
    fw3 = gen_relax_fw(
        atoms,
        db_name_remote,
        workdir + "/tight_relax/",
        "light_relax_atoms",
        "tight_relax_atoms",
        db_label="tight_relax",
        up_calc_from_db=["k_grid"],
        name=f"tight_relax_{name}",
        spec_qad=spec_qad_relax,
        from_db=True,
        calc_settings=relax_tight,
        hilde_cfg=hilde_cfg
    )
    fw4 = gen_initialize_phonopy_fw(
        atoms,
        smatrix,
        workdir + "/force_calcs/",
        "tight_relax_atoms",
        symprec=symprec,
        up_calc_from_db=["k_grid"],
        name=f"init_phono_{name}",
        spec_qad=spec_qad_forces,
        from_db=True,
        calc_settings=force_calc,
        hilde_cfg=hilde_cfg
    )
    fw5 = gen_analyze_phonopy_fw(
        atoms,
        db_name_local,
        smatrix,
        "tight_relax_atoms",
        db_label="phonons",
        symprec=symprec,
        name=f"analyze_phono_{name}",
        from_db=True
    )
    workflow = Workflow([fw1, fw2, fw3, fw4, fw5], {fw1:[fw2], fw2:[fw3], fw3:[fw4], fw4:[fw5]},
                        name=name)
    return workflow
