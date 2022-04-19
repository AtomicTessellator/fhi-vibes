def cond_test_atomic_volume(atoms):
    """Returns true if atomic volume is less than 100.0 AA^3"""
    return atoms.get_volume() / len(atoms) < 100.0
