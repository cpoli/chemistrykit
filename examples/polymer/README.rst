Examples
========

This gallery walks through every public feature of
``chemistrykit.polymer``: ideal random-walk chain statistics and the
Flory exponent for real chains under different solvent conditions;
molecular-weight-distribution statistics and the closed-form Flory-Schulz
distribution; step-growth kinetics via the Carothers equation; and
chain-growth/free-radical polymerization kinetics built on
``chemistrykit.kinetics``'s reaction-network engine.

Each script in this gallery is self-contained and can be run directly
with ``python examples/polymer/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **chain_statistics** -- ideal-chain end-to-end distance and radius of
  gyration, and the Flory-exponent scaling of real chains under
  theta/good/poor solvent conditions.
- **molecular_weight_distribution** -- Mn/Mw/PDI, and the closed-form
  Flory-Schulz (most-probable) chain-length distribution, cross-checked
  by direct numerical summation.
- **step_growth** -- the Carothers equation relating degree of
  polymerization to extent of reaction, plain and with a stoichiometric
  imbalance.
- **chain_growth** -- free-radical initiation/propagation/termination
  kinetics, integrated numerically and checked against the steady-state
  approximation's closed form.
