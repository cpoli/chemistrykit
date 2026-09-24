Examples
========

This gallery walks through every public feature of ``chemistrykit.md``:
the Lennard-Jones fluid in reduced units (energy conservation, pressure,
and the radial distribution function g(r)); Morse/Buckingham/harmonic
bonded potentials and small bonded clusters; periodic boundaries; SHAKE
bond constraints; thermostats; and self-diffusion coefficients.

Each script in this gallery is self-contained and can be run directly with
``python examples/md/<section>/<script>.py``. Every script also carries an
RST module docstring as its title/description and uses ``# %%`` markers to
split narrative text from code, which is exactly what Sphinx-Gallery
renders into the pages below -- the script *is* the source of truth for
what you see, not a copy of it.

Sections
--------

- **lj_fluid** -- Lennard-Jones fluid simulations: Rahman's liquid argon
  g(r), a cross-check of the simulated speed distribution against
  :mod:`chemistrykit.statmech`'s Maxwell-Boltzmann distribution,
  hard-sphere-like packing order, the Weeks-Chandler-Andersen repulsive
  reference fluid, Verlet's integrator and neighbor list, and Einstein
  and Green-Kubo self-diffusion coefficients.
- **pair_potentials** -- the Lennard-Jones potential, Morse vs. harmonic
  bond potentials and the classical vibration of a two-body
  :class:`~chemistrykit.md.systems.pair_potentials.DiatomicOscillator`,
  the Buckingham potential, and a small bonded molecule vibrating under
  coupled harmonic bond-stretch and angle-bend terms.
- **periodic_boundaries** -- periodic boundary conditions and the
  minimum-image convention.
- **constraints** -- rigid bonds enforced with the SHAKE algorithm.
- **thermostats** -- Nose-Hoover, Berendsen weak-coupling, and stochastic
  velocity-rescaling temperature control.
