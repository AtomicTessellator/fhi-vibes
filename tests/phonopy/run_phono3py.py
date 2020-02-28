""" use the vibes phonopy workflow """


def main():
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from vibes.phono3py import run_phono3py

    atoms = bulk("Al")

    calc = EMT()

    run_phono3py(atoms=atoms, calculator=calc)


try:
    import phono3py  # noqa: F401

    main()

except ModuleNotFoundError:
    print("phono3py not installed, skip")
