Examples
========

This gallery walks through every public feature of
``chemistrykit.photochem``: Jablonski-diagram excited-state kinetics
(reusing ``chemistrykit.kinetics``'s reaction-network engine),
fluorescence/phosphorescence quantum yields and the photochemical
quantum yield via Beer-Lambert, Stern-Volmer quenching with a
static-vs-dynamic diagnostic, and photostationary-state kinetics for a
two-state photoswitch.

Each script in this gallery is self-contained and can be run directly
with ``python examples/photochem/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **jablonski** -- the 3-state Jablonski excited-state decay network,
  checked against its closed-form population solution.
- **quantum_yield** -- fluorescence/phosphorescence quantum yields, and
  the photochemical quantum yield via Beer-Lambert photon absorption.
- **stern_volmer** -- Stern-Volmer quenching, fitting a quenching
  constant, and distinguishing static from dynamic quenching.
- **photostationary_state** -- the photostationary state of a two-state
  photoswitch under simultaneous forward/reverse photolysis, checked
  against long-time numerical integration.
