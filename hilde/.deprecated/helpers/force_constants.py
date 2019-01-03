""" Helpers to convert force constants to dynamical matrix etc """


def fc2dynmat(atoms, force_constants):
    """ Convert force_constants for atoms into dynamical matrix """

    masses = atoms.get_masses()
    rminv = (masses ** -0.5).repeat(3)

    dynamical_matrix = force_constants * rminv[:, None] * rminv[None, :]

    return dynamical_matrix

