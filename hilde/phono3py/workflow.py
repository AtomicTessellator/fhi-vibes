""" Provide a highlevel phono3py workflow for computing 3rd order force constants """

from hilde.tasks import calculate_socket
from hilde.helpers.restarts import restart
from hilde.phonopy.workflow import bootstrap

from .postprocess import postprocess


def run_phono3py(postprocess_args=None, **kwargs):
    """ high level function to run phono3py workflow

    Parameters
    ----------
    postprocess_args: dict
        arguments for postprocessing
    kwargs: dict
        Preprocessing arguments, items must include

        atoms: ase.atoms.Atoms
            primitive cell for the calculation
        supercell_matrix: np.ndarray
            supercell matrix for the third order phonons

        items may include

        cutoff_pair_distance: float
            All pairs further apart than this cutoff are ignored
        is_diagonal: bool
            Whether allow diagonal displacements of Atom 2 or not
        q_mesh: np.ndarray
            q-point interpolation mesh postprocessing
        displacement: float
            magnitude of the displacement
        symprec: float
            distance tolerance for determining the sapce group/symmetry
        log_level: int
            How much information should be streamed to the console
    """

    args = bootstrap(name="phono3py", **kwargs)

    completed = calculate_socket(**args)

    if not completed:
        restart()
    else:
        print("Start postprocess.")
        if postprocess_args is None:
            postprocess_args = {}
        args.update(postprocess_args)
        postprocess(**args)
        print("done.")
