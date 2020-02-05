"""Naming the dimensions

    Atom labels: I, J
    Cartesian coordinates: a, b
"""
I, J, a, b = "I", "J", "a", "b"

# composite
Ia, Jb = "Ia", "Jb"

time = "time"
atoms = [time, I]
vec = [time, a]
atoms_vec = [time, I, a]
stress = [time, a, b]
stresses = [time, I, a, b]
kappa = [time, I, J, a, b]

lattice = [a, b]
positions = [I, a]

fc_remapped = [Ia, Jb]
