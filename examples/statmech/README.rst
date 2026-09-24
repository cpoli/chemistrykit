Examples
========

This gallery walks through every public feature of
``chemistrykit.statmech``: translational/rotational/vibrational partition
functions and the thermodynamic functions derived from them; the
Maxwell-Boltzmann speed distribution; the lattice-gas adsorption model;
quantum statistics and its classical limit; the Debye solid; exact
Ising-model results; and the second virial coefficient.

Each script in this gallery is self-contained and can be run directly with
``python examples/statmech/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **partition_functions** -- equipartition, the Einstein vibrational
  heat capacity, Gibbs's canonical partition function, and the
  Sackur-Tetrode translational entropy.
- **maxwell_boltzmann** -- the Maxwell-Boltzmann speed distribution, its
  characteristic speeds, and a Stern-type molecular-beam check.
- **lattice_gas** -- the lattice-gas derivation of the Langmuir adsorption
  isotherm, and :math:`S=k_B\ln W` by counting arrangements.
- **quantum_statistics** -- the thermal de Broglie wavelength, and the
  Maxwell-Boltzmann limit of Bose-Einstein and Fermi-Dirac statistics.
- **solids** -- Debye's :math:`T^3` law for the heat capacity of a solid.
- **ising** -- Ising's 1D chain, Kramers-Wannier duality, and Onsager's
  exact 2D solution.
- **virial** -- Mayer's second virial coefficient from a pair potential.
