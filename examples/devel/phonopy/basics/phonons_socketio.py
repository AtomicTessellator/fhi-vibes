from time import time
from pathlib import Path
import matplotlib.pyplot as plt

from vibes.parsers import read_structure
from vibes.helpers.supercell import make_cubic_supercell
from vibes.helpers import d2k
from vibes.phonopy import wrapper as ph
from vibes.tasks.calculate import calculate_multiple_socketio
from vibes.templates.aims import setup_aims
import vibes.helpers.brillouinzone as bz

atoms = read_structure("../si.in")
vol = atoms.get_volume()

base_folder = Path(f"{atoms.sysname}").absolute()

port = 31415
k_grid = d2k(atoms, 1)
n_target = 4

print(f"k_grid primitive cell: {k_grid}")

relax_folder = base_folder / f"relax_{vol:.3f}"

try:
    ratoms = read_structure(relax_folder / "geometry.in.next_step")
except FileNotFoundError:
    ratoms = atoms

sc, smatrix = make_cubic_supercell(ratoms, target_size=n_target)

print(sc.cell)

phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

workdir = base_folder / "{}{}{}_{}{}{}_{}{}{}_{:.3f}".format(*smatrix.flatten(), vol)
workdir.mkdir(exist_ok=True, parents=True)

sc.write(workdir / "supercell.in")

k_grid = d2k(sc, 1)
print(f"k_grid supercell: {k_grid}")

# set the phonon object
phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

print(f"{len(scs)} supercells created.")

socketio_settings = {
    "use_pimd_wrapper": ("localhost", port),
    "compute_forces": True,
    "sc_accuracy_rho": 1e-4,
    "k_grid": k_grid,
}

aims = setup_aims(custom_settings=socketio_settings)

phonopy_log = workdir / "phonopy.log"
trajectory = workdir / "phonopy.traj"
tmp_dir = workdir / "socketio"

stime = time()
cells_computed = calculate_multiple_socketio(
    cells=scs,
    calculator=aims,
    workdir=tmp_dir,
    port=port,
    trajectory=trajectory,
    force=False,
    log_file=phonopy_log,
)
timing_socket_io = time() - stime
print(f".. done in {timing_socket_io:.2f}s")

force_sets = [atoms.get_forces() for atoms in cells_computed]
phonon.produce_force_constants(force_sets)

dos = ph.get_dos(phonon)
plt.plot(dos[0], dos[1])
plt.savefig(workdir / "dos.pdf")

# bands
bands, labels = bz.get_bands_and_labels(phonon.primitive)
phonon.set_band_structure(bands)
phonon.plot_band_structure(labels=labels)
plt.savefig(workdir / "bandstructure.pdf")
