"""Naming the dimensions

    Atom labels: I, J
    Cartesian coordinates: a, b
"""
I, J, a, b = "I", "J", "a", "b"

# composite
Ia, Jb = "Ia", "Jb"

time = "time"
time_atom = (time, I)
time_vec = (time, a)
time_atom_vec = (time, I, a)
time_tensor = (time, a, b)
# time_atom_tensor = (time, I, a, b)
# time_atom_atom_tensor = (time, I, J, a, b)
# taat = time_atom_atom_tensor

lattice = (a, b)
positions = (I, a)

fc_remapped = (Ia, Jb)
