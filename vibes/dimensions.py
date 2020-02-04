"""Naming the dimensions"""
I, J, a, b = "I", "J", "a", "b"

time = "time"
atoms = [time, I]
vec = [time, a]
atoms_vec = [time, I, a]
stress = [time, a, b]
stresses = [time, I, a, b]
kappa = [time, I, J, a, b]

lattice = [a, b]
positions = [I, a]
