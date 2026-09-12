Examples
========

This gallery walks through every public feature of
``chemistrykit.surface``: Langmuir, Freundlich, and BET adsorption
isotherms with their standard linearizations for fitting parameters from
data; Langmuir-Hinshelwood single- and dual-site surface-reaction
kinetics; and a turnover-frequency/rate-enhancement catalysis model
built on ``chemistrykit.kinetics``'s Arrhenius equation.

Each script in this gallery is self-contained and can be run directly
with ``python examples/surface/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **langmuir** -- the Langmuir monolayer isotherm, its exact
  half-saturation point, and its linearization for fitting.
- **freundlich** -- the empirical Freundlich power-law isotherm and its
  linearization.
- **bet** -- the BET multilayer isotherm, and its numerical reduction to
  Langmuir in the appropriate limit.
- **langmuir_hinshelwood** -- single- and dual-site surface-reaction
  kinetics built on the Langmuir coverage.
- **catalysis** -- turnover frequency/number, and rate enhancement from a
  catalyst's activation-energy reduction.
