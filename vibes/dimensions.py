"""Naming the dimensions

    Atom labels: I, J
    Cartesian coordinates: a, b
    Mode labels: s, q
"""
i, j, I, J, a, b = "i", "j", "I", "J", "a", "b"
s, q = "s", "q"
s_q = (s, q)
q_a = (q, a)

# composite
Ia, Jb = "Ia", "Jb"
ia, jb = "ia", "jb"

vec = (a,)
tensor = (a, b)
ab = tensor

time = "time"
time_atom = (time, I)
time_vec = (time, *vec)
time_atom_vec = (time, I, a)
time_tensor = (time, *tensor)
time_atom_tensor = (time, I, *tensor)
# time_atom_atom_tensor = (time, I, J, a, b)
# taat = time_atom_atom_tensor

lattice = tensor
positions = (I, a)

fc = (i, J, a, b)
fc_remapped = (Ia, Jb)

dim_replace = {"a1": a, "a2": b, "I1": I, "I2": J}
