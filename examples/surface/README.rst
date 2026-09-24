Examples
========

This gallery walks through every public feature of
``chemistrykit.surface``: Langmuir, Freundlich, and BET adsorption
isotherms with their standard linearizations for fitting parameters from
data; Langmuir-Hinshelwood single- and dual-site surface-reaction
kinetics; and a turnover-frequency/rate-enhancement catalysis model
built on ``chemistrykit.kinetics``'s Arrhenius equation; the Gibbs
adsorption equation; Polanyi/Dubinin-Radushkevich and Temkin isotherms;
Eley-Rideal kinetics; and temperature-programmed desorption.

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
- **catalysis** -- one example per catalysis milestone: Doebereiner and
  Berzelius, Arrhenius and Ostwald, Sabatier's hydrogenation, Constable's
  compensation effect, Balandin's volcano curve, and Boudart's turnover
  frequency.
- **gibbs_adsorption** -- surface excess from surface tension via the
  Gibbs adsorption equation.
- **dubinin** -- Polanyi's adsorption potential and the
  Dubinin-Radushkevich characteristic curve.
- **temkin** -- the logarithmic Temkin isotherm from a spread of site
  energies.
- **eley_rideal** -- Eley-Rideal kinetics contrasted with dual-site
  Langmuir-Hinshelwood kinetics.
- **tpd** -- temperature-programmed desorption and Redhead's analysis.
