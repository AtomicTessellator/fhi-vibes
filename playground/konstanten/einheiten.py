from numpy import pi
Bohr_to_AA       =  0.52917721092
Hartree_to_eV    = 27.211385
#u_in_kg          =  1.6726219e-27
u_in_kg          =  1.66053904e-27
au_to_kg         =  9.10938291e-31
kg_to_au         =  1 / au_to_kg
u_in_au          =  u_in_kg * kg_to_au
c_in_au          =  137
au_to_s          =  2.418884326505e-17
au_to_fs         =  2.418884326505e-2
au_to_cm         =  au_to_s * 2*pi * 3e10
eV_to_THz        = 15.633302  #241.8 =98.3/2/pi
au_to_K          = 3.1577464e5
THz_to_cm        = 33.36
eV_to_cm         = eV_to_THz * THz_to_cm

# physical constants
EV = 1.60217733e-19 # [J]
Avogadro = 6.02214179e23
PlanckConstant = 4.13566733e-15 # [eV s]
kb_J = 1.3806504e-23 # [J/K]

Kb = kb_J / EV  # [eV/K] 8.6173383e-05
EvTokJmol = EV / 1000 * Avogadro # [kJ/mol] 96.4853910
kJmolToEv = 1 / EvTokJmol
THzToEv = PlanckConstant * 1e12 # [eV]