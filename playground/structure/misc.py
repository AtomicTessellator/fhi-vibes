def generate_lattice(a, b=None, c=None, alpha=90, beta=90, gamma=90, lattice_type=None):
    """ [pymatgen, adapted]
    Create a Lattice using unit cell lengths (Angstrom) and angles (in degrees).

    Args:
        cellpars (list of floats):
            a (float): *a* lattice parameter.
            b (float): *b* lattice parameter.
            c (float): *c* lattice parameter.
            alpha (float): *alpha* angle in degrees.
            beta (float): *beta* angle in degrees.
            gamma (float): *gamma* angle in degrees.
        lattice_type (str):

    Returns:
        Lattice cell matrix with the specified lattice parameters.
    """
    if lattice_type == 'cubic':
        return np.array( [[a, 0.0, 0.0],
                          [0.0, a, 0.0],
                          [0.0, 0.0, a]] )
    elif lattice_type == 'tetragonal':
        if c is None:
            print( 'Error: Tetragonal lattice needs parameter `c`.' )
            return None
        b = a
    elif lattice_type == 'orthorhombic':
        if b is None or c is None:
            print( 'Error: Orthorhombic lattice needs parameters `b` and `c`.' )
    elif lattice_type == 'monoclinic':
        if b is None or c is None:
            print( 'Error: Monoclinic lattice needs parameters `b` and `c`.' )
        if beta == 90:
            print( 'Warning: Have you set beta? It is 90 -> orthorhombic.' )
    elif lattice_type == 'hexagonal':
        if c is None:
            print( 'Error: Hexagonal lattice needs parameters `c`.' )
        b = a
        gamma = 120
    elif lattice_type == 'rhombohedral':
        b = a
        c = a
        beta = alpha
        gamma = alpha

    alpha_r = np.radians(alpha)
    beta_r = np.radians(beta)
    gamma_r = np.radians(gamma)
    val = (cos(alpha_r) * cos(beta_r) - cos(gamma_r)) \
          / (sin(alpha_r) * sin(beta_r))
    # Sometimes rounding errors result in values slightly > 1.
    val = max(min(val, 1), -1)
    gamma_star = np.arccos(val)
    vector_a = [a * sin(beta_r), 0.0, a * cos(beta_r)]
    vector_b = [-b * sin(alpha_r) * cos(gamma_star),
                b * sin(alpha_r) * sin(gamma_star),
                b * cos(alpha_r)]
    vector_c = [0.0, 0.0, float(c)]
    return np.array([vector_a, vector_b, vector_c])
