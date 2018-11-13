aims_kgrid_conv_settings = {
    "species_type" : "light",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_iter_limit": 1000
}
aims_relax_settings_light = {
    "species_type" : "light",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "relax_geometry": "trm 1E-2",
    "relax_unit_cell": "full",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}
aims_relax_settings_tight = {
    "species_type" : "tight",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "relax_geometry": "trm 5E-3",
    "relax_unit_cell": "full",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}
aims_force_settings = {
    "species_type" : "tight",
    "output_level": 'MD_light',
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "include_spin_orbit": "non_self_consistent",
    "symmetry_reduced_k_grid": "false.",
    "occupation_type": "gaussian 0.01",
    "mixer": "pulay",
    "n_max_pulay": 8,
    "charge_mix_param": 0.5,
    "sc_accuracy_eev": 1E-3,
    "sc_accuracy_rho": 1E-5,
    "sc_accuracy_etot": 1E-6,
    "sc_accuracy_forces": 5E-4,
    "sc_iter_limit": 1000
}