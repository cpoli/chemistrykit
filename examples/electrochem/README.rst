Examples
========

This gallery walks through every public feature of
``chemistrykit.electrochem``: the Nernst equation for standard and
concentration cells (with activity-coefficient corrections), a curated
standard-reduction-potential table with redox-couple balancing,
Butler-Volmer electrode kinetics and Tafel-plot linearization, Faraday's
laws of electrolysis and the galvanic-vs-electrolytic distinction, a
simplified constant-current battery discharge model with Peukert's-law
rate dependence, fuel-cell thermodynamics, Kohlrausch's conductivity
laws, and diffusion-limited electroanalytical currents.

Each script in this gallery is self-contained and can be run directly
with ``python examples/electrochem/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **nernst** -- the Nernst equation for standard cells and concentration
  cells, plus Debye-Huckel activity-coefficient-corrected reaction
  quotients.
- **standard_potentials** -- the curated standard-reduction-potential
  table, redox-couple electron balancing, and cell-potential combination.
- **butler_volmer** -- Butler-Volmer electrode kinetics and its
  high-overpotential Tafel-plot linearization, checked for convergence.
- **electrolysis** -- Faraday's laws of electrolysis, and galvanic vs.
  electrolytic cell operation.
- **battery** -- Volta's pile, a simplified constant-current battery
  discharge model, and Peukert's-law capacity-vs-rate dependence.
- **fuel_cell** -- Grove's hydrogen-oxygen gas battery and its
  thermodynamic voltage and efficiency limits.
- **conductivity** -- Kohlrausch's laws of electrolytic conductivity.
- **voltammetry** -- the Cottrell equation, polarography and the Ilkovič
  equation, and the Randles-Ševčík peak current.
