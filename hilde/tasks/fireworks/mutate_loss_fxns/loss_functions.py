module_name = __name__

def energy_diff( atoms_cur, atoms_prev, criteria=0.01 ):
    ener_diff = atoms_cur["results"]["energy"] - atoms_prev["results"]["energy"]
    return abs(ener_diff / atoms_cur["results"]["energy"]) < criteria

energy_diff.name = f'{module_name}.{energy_diff.__name__}'
