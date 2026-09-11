Examples
========

This gallery walks through every public feature of ``chemistrykit.md``:
the Lennard-Jones fluid in reduced units (energy conservation, pressure,
and the radial distribution function g(r)); Morse/Buckingham/harmonic
bonded potentials and small bonded clusters; and velocity-rescaling and
Nose-Hoover thermostats.

Each script in this gallery is self-contained and can be run directly with
``python examples/md/<section>/<script>.py``. Every script also carries an
RST module docstring as its title/description and uses ``# %%`` markers to
split narrative text from code, which is exactly what Sphinx-Gallery
renders into the pages below -- the script *is* the source of truth for
what you see, not a copy of it.

Sections
--------

- **lj_fluid** -- an NVE Lennard-Jones fluid run: energy conservation,
  the radial distribution function g(r), and a cross-check of the
  simulated speed distribution against
  :mod:`chemistrykit.statmech`'s Maxwell-Boltzmann distribution.
- **pair_potentials** -- Morse vs. harmonic bond potentials, and the
  classical vibration of a two-body :class:`~chemistrykit.md.systems.pair_potentials.DiatomicOscillator`.
- **thermostats** -- comparing an unthermostatted (NVE) run against
  velocity-rescaling and Nose-Hoover temperature control.
