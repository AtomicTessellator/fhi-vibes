""" Provide a full highlevel phonopy workflow """
from pathlib import Path

from phonopy.file_IO import write_FORCE_CONSTANTS

from hilde.helpers.brillouinzone import get_special_points
from hilde.helpers.converters import dict2atoms
from hilde.helpers import Timer
from hilde.helpers.paths import cwd
from hilde.phonopy.wrapper import (
    prepare_phonopy,
    get_force_constants,
    plot_bandstructure as plot_bs,
    get_bandstructure,
    plot_bandstructure_and_dos,
    get_animation,
)
from hilde.phonopy import defaults
from hilde.structure.convert import to_Atoms, to_Atoms_db
from hilde.trajectory import reader
from hilde.io import write
from hilde.helpers import warn, talk, Timer

from . import displacement_id_str


def postprocess(
    trajectory="trajectory.son",
    workdir=".",
    calculate_full_force_constants=False,
    born_charges_file=None,
    verbose=True,
    **kwargs,
):
    """ Phonopy postprocess """

    timer = Timer()

    trajectory = Path(workdir) / trajectory

    talk("Start phonopy postprocess:")

    calculated_atoms, metadata = reader(trajectory, True)

    # make sure the calculated atoms are in order
    for nn, atoms in enumerate(calculated_atoms):
        atoms_id = atoms.info[displacement_id_str]
        if atoms_id == nn:
            continue
        warn(f"Displacement ids are not in order. Inspect {trajectory}!", level=2)

    for disp in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
        disp["number"] = int(disp["number"])
    primitive = dict2atoms(metadata["Phonopy"]["primitive"])
    supercell = dict2atoms(metadata["atoms"])
    supercell_matrix = metadata["Phonopy"]["supercell_matrix"]
    supercell.info = {"supercell_matrix": str(supercell_matrix)}
    symprec = metadata["Phonopy"]["symprec"]

    phonon = prepare_phonopy(primitive, supercell_matrix, symprec=symprec)
    phonon._displacement_dataset = metadata["Phonopy"]["displacement_dataset"].copy()

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    phonon.produce_force_constants(
        force_sets, calculate_full_force_constants=calculate_full_force_constants
    )

    # born charges?
    if born_charges_file:
        from phonopy.file_IO import get_born_parameters

        prim = phonon.get_primitive()
        psym = phonon.get_primitive_symmetry()
        if verbose:
            talk(f".. read born effective charges from {born_charges_file}")
        nac_params = get_born_parameters(open(born_charges_file), prim, psym)
        phonon.set_nac_params(nac_params)

    if calculate_full_force_constants:
        phonon.produce_force_constants(force_sets, calculate_full_force_constants=True)

        # force_constants = get_force_constants(phonon)
        # fname = "force_constants.dat"
        # np.savetxt(fname, force_constants)
        # talk(f".. Force constants saved to {fname}.")

    if verbose:
        timer("done")
    return phonon


def extract_results(
    phonon,
    write_geometries=True,
    write_force_constants=True,
    write_thermal_properties=False,
    write_bandstructure=False,
    write_dos=False,
    write_pdos=False,
    plot_bandstructure=True,
    plot_dos=False,
    plot_pdos=False,
    animate=None,
    animate_q=None,
    q_mesh=None,
    output_dir="phonopy_output",
    tdep=False,
    tdep_reduce_fc=True,
):
    """ Extract results from phonopy object and present them.
        With `tdep=True`, the necessary input files for TDEP's
          `convert_phonopy_to_forceconstant`
        are written. """

    timer = Timer("\nExtract phonopy results:")
    if q_mesh is None:
        q_mesh = defaults.q_mesh.copy()
    talk(f".. q_mesh:   {q_mesh}")

    primitive = to_Atoms(phonon.get_primitive())
    supercell = to_Atoms(phonon.get_supercell())

    Path.mkdir(Path(output_dir), exist_ok=True)
    with cwd(output_dir):
        if write_geometries:
            talk(f".. write primitive cell")
            write(primitive, "geometry.in.primitive")
            talk(f".. write supercell")
            write(supercell, "geometry.in.supercell")

        if write_force_constants:
            talk(f".. write force constants")
            write_FORCE_CONSTANTS(
                phonon.get_force_constants(),
                filename="FORCE_CONSTANTS",
                p2s_map=phonon.get_primitive().get_primitive_to_supercell_map(),
            )

        if write_thermal_properties:
            talk(f".. write thermal properties")
            phonon.run_mesh(q_mesh)
            phonon.run_thermal_properties()
            phonon.write_yaml_thermal_properties()

        if plot_bandstructure:
            talk(f".. plot band structure")
            plot_bs(phonon, file="bandstructure.pdf")
        if write_bandstructure:
            talk(f".. write band structure yaml file")
            get_bandstructure(phonon)
            phonon.write_yaml_band_structure()

        if write_dos:
            talk(f".. write DOS")
            phonon.run_mesh(q_mesh, with_eigenvectors=True)
            phonon.run_total_dos(use_tetrahedron_method=True)
            phonon.write_total_dos()
        if plot_dos:
            talk(f".. plot DOS")
            plot_bandstructure_and_dos(phonon, file="bands_and_dos.pdf")

        if write_pdos:
            talk(f".. write projected DOS")
            phonon.run_mesh(q_mesh, with_eigenvectors=True, is_mesh_symmetry=False)
            phonon.run_projected_dos(use_tetrahedron_method=True)
            phonon.write_projected_dos()

        if plot_pdos:
            talk(f".. plot projected DOS")
            plot_bandstructure_and_dos(phonon, partial=True, file="bands_and_pdos.pdf")

        animate_q_points = []
        if animate:
            for q_pt in get_special_points(primitive).values():
                animate_q_points.append(tuple(q_pt))

        elif animate_q:
            animate_q_points = animate_q

        for q_pt in animate_q_points:
            path = Path("animation")
            path.mkdir(exist_ok=True)
            outfile = path / f"animation_{q_pt[0]}_{q_pt[1]}_{q_pt[2]}.ascii"
            get_animation(phonon, q_pt, outfile)
            talk(f".. {outfile} written")

    if tdep:
        write_settings = {"format": "vasp", "direct": True, "vasp5": True}
        fnames = {"primitive": "infile.ucposcar", "supercell": "infile.ssposcar"}

        # reproduce reduces force constants
        if tdep_reduce_fc:
            phonon.produce_force_constants(calculate_full_force_constants=False)

            fname = fnames["primitive"]
            primitive.write(fname, **write_settings)
            talk(f"Primitive cell written to {fname}")

            fname = fnames["supercell"]
            supercell.write(fname, **write_settings)
            talk(f"Supercell cell written to {fname}")

    timer(f"all files written to {output_dir}")
