Documentation
===

Here we give an overview of the features of `FHI-vibes` and document the basic usage principles.

In particular, you will find

- [Units](units.md)
- [Input files](input_files.md)
- [Output files](output_files.md)
- [Calculator setup](calculator_setup.md)
- Tasks and workflows, i.e.,
    - [Geometry optimization](relaxation.md)
    - [Phonon calculations](phonopy.md)
    - [Phono3py calculations](phono3py.md)
    - [Molecular dynamics simulations](md.md)
- and [High-Throughput workflows.](../High_Throughput/Documentation/0_overview.md)

Detailed introductions into specific tasks can be found in the [Tutorial](../Tutorial/0_intro.md).



## Design Philosophy

The general design principles are

- calculations should be defined with human-readable input files,
- a calculation should produce a self-contained output file that _includes_ metadata describing the calculation, as well as the calculated results.
- Output files should be easy to parse by a computer.
