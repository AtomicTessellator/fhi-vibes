default_settings = {
    "species_type": "light",
    "output_level": "MD_light",
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "sc_accuracy_rho": 1e-5,
}

aims_kgrid_conv_settings = default_settings

aims_relax_settings_light = {
    **default_settings,
    "relax_geometry": "lattice_trm 1E-2",
    "relax_unit_cell": "full",
    "sc_accuracy_rho": 1e-6,
}

aims_relax_settings_tight = {
    **aims_relax_settings_light,
    "species_type": "tight",
    "relax_geometry": "lattice_trm 1E-3",
}

aims_force_settings = {
    **default_settings,
    "species_type": "tight",
    "sc_accuracy_rho": 1e-6,
    "compute_forces": True,
}
