""" Provide a full highlevel phonopy workflow """
from pathlib import Path

from phonopy.file_IO import write_FORCE_CONSTANTS

from hilde.helpers.brillouinzone import get_special_points
from hilde.helpers.converters import dict2atoms
from hilde.helpers.paths import cwd
from hilde.phonopy import wrapper
from hilde.phonopy import defaults
from hilde.structure.convert import to_Atoms
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
    """Phonopy postprocess

    Parameters
    ----------
    trajectory: str or Path
        The trajectory file to process
    workdir: str or Path
        The working directory where trajectory is stored
    calculate_full_force_constants: bool
        If True calculate the full force constant matrix
    born_charges_file: str or Path
        Path to the born charges file
    verbose: bool
        If True be verbose

    Returns
    -------
    phonon: phonopy.Phonopy
        The Phonopy object with the force constants calculated
    """

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

    phonon = wrapper.prepare_phonopy(primitive, supercell_matrix, symprec=symprec)
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
    plot_thermal_properties=False,
    plot_dos=False,
    plot_pdos=False,
    animate=None,
    animate_q=None,
    q_mesh=None,
    output_dir="phonopy_output",
    tdep=False,
    tdep_reduce_fc=True,
    verbose=False,
):
    """ Extract results from phonopy object and present them.

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The Phonopy Object with calculated force constants
    write_geometries: bool
        If True write the geometry files for the primitive and supercells
    write_force_constants: bool
        If True write the FORCE_CONSTANTS file
    write_thermal_properties: bool
        If True write a thermal properties yaml file
    write_bandstructure: bool
        If True write the band.yaml file
    write_dos: bool
        If True write the total DOS output file
    write_pdos: bool
        If True the projected DOS output file
    plot_bandstructure: bool
        If True plot the band structure save it to a pdf
    plot_thermal_properties: bool
        If True plot the thermal properties and save them to a pdf
    plot_dos: bool
        If True plot the total density of states and save it to a pdf
    plot_pdos: bool
        If True plot the projected density of states and save it to a pdf
    animate: bool
        If True write anaimation files for all high-symmetry poitns
    animate_q: list of tuples
        A list of q_points to write animation files for
    q_mesh: np.ndarray
        The size of the interpolated q-grid
    output_dir: str or Path
        Directory to store output files
    tdep: bool
        If True the necessary input files for TDEP's `convert_phonopy_to_forceconstant` are written.
    tdep_reduce_fc: bool
        If True reduce force_constants to match tdep's format
    """

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
        if plot_thermal_properties:
            talk(f".. plot thermal properties")
            wrapper.plot_thermal_properties(phonon)

        if plot_bandstructure:
            talk(f".. plot band structure")
            wrapper.plot_bandstructure(phonon, file="bandstructure.pdf")
        if write_bandstructure:
            talk(f".. write band structure yaml file")
            wrapper.set_bandstructure(phonon)
            phonon.write_yaml_band_structure()

        if write_dos:
            talk(f".. write DOS")
            phonon.run_mesh(q_mesh, with_eigenvectors=True)
            phonon.run_total_dos(use_tetrahedron_method=True)
            phonon.write_total_dos()
        if plot_dos:
            talk(f".. plot DOS")
            wrapper.plot_bandstructure_and_dos(phonon, file="bands_and_dos.pdf")

        if write_pdos:
            talk(f".. write projected DOS")
            phonon.run_mesh(q_mesh, with_eigenvectors=True, is_mesh_symmetry=False)
            phonon.run_projected_dos(use_tetrahedron_method=True)
            phonon.write_projected_dos()

        if plot_pdos:
            talk(f".. plot projected DOS")
            wrapper.plot_bandstructure_and_dos(
                phonon, partial=True, file="bands_and_pdos.pdf"
            )

        animate_q_points = {}
        if animate:
            animate_q_points = get_special_points(primitive)

        elif animate_q:
            for q_pt in animate_q:
                key = "_".join(str(q) for q in q_pt)
                animate_q_points.update({f"{key}": q_pt})

        for key, val in animate_q_points.items():
            path = Path("animation")
            path.mkdir(exist_ok=True)
            outfile = path / f"animation_{key}.ascii"
            wrapper.get_animation(phonon, val, outfile)
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

    if verbose:
        talk("\nFrequencies at Gamma point:")
        phonon.run_mesh([1, 1, 1])
        qpoints, weights, frequencies, _ = phonon.get_mesh()
        for q, w, f in zip(qpoints, weights, frequencies):
            print(f"q = {q} (weight= {w})")
            print("# Mode   Frequency")
            for ii, fi in enumerate(f):
                print(f"  {ii+1:3d} {fi:12.7f} THz")
