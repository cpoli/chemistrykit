Examples
========

This gallery walks through every public feature of
``chemistrykit.photochem``: Jablonski-diagram excited-state kinetics
(reusing ``chemistrykit.kinetics``'s reaction-network engine),
fluorescence/phosphorescence quantum yields and the photochemical
quantum yield via Beer-Lambert, Stern-Volmer quenching with a
static-vs-dynamic diagnostic, photostationary-state kinetics for a
two-state photoswitch, the H2/Cl2 photochemical chain reaction,
fluorescence observables (Stokes shift, Perrin anisotropy), Förster and
Dexter energy transfer, ferrioxalate actinometry, and Rehm-Weller
electron-transfer quenching.

Each script in this gallery is self-contained and can be run directly
with ``python examples/photochem/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **jablonski** -- the 3-state Jablonski excited-state decay network,
  checked against its closed-form population solution; triplet
  phosphorescence, flash photolysis of a transient triplet, and Kasha's
  rule.
- **quantum_yield** -- the Grotthuss-Draper and Stark-Einstein laws of
  photochemistry, and Vavilov's excitation-independent fluorescence
  quantum yield.
- **chain_reaction** -- the H2/Cl2 photochemical chain and its
  super-unity quantum yield.
- **stern_volmer** -- Stern-Volmer quenching, fitting a quenching
  constant, and distinguishing static from dynamic quenching.
- **photostationary_state** -- the photostationary state of a two-state
  photoswitch under simultaneous forward/reverse photolysis, checked
  against long-time numerical integration.
- **fluorescence** -- the Stokes shift and the Perrin anisotropy equation.
- **energy_transfer** -- Förster (FRET) and Dexter exchange energy transfer.
- **actinometry** -- photon-flux measurement with the ferrioxalate
  actinometer.
- **electron_transfer** -- Rehm-Weller electron-transfer quenching.
