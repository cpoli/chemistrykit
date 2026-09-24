Examples
========

This gallery walks through every public feature of
``chemistrykit.crystal``: the 7 crystal systems and general unit-cell
volume; hard-sphere packing (packing fraction, coordination number) for
the SC/BCC/FCC/HCP lattices; ionic-crystal lattice energy via the
Born-Lande and Kapustinskii equations, backed by a genuinely converging
(Evjen-method) numerical Madelung constant; Bragg's law and powder-XRD
peak positions with structure factors and systematic absences, Laue
interference, and Scherrer crystallite sizes; Schottky/Frenkel
point-defect equilibrium; Miller indices, interfacial angles, and the 14
Bravais lattices; and Pauling's radius-ratio rule and Goldschmidt's
perovskite tolerance factor.

Each script in this gallery is self-contained and can be run directly
with ``python examples/crystal/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **crystal_systems** -- classifying a unit cell into one of the 7
  crystal systems, the general lattice-parameter volume formula, Steno's
  interfacial angles, Miller indices from face intercepts, and the 14
  Bravais lattices.
- **packing** -- atomic packing factor, coordination number, and atoms
  per cell for SC, BCC, FCC, and ideal HCP.
- **lattice_energy** -- the Born-Lande and Kapustinskii equations, and
  Shannon's effective ionic radii.
- **madelung** -- Madelung's lattice sum shell by shell, and the NaCl
  Madelung constant from a genuinely converging
  (Evjen-method) lattice summation, contrasted with a naive truncated sum
  that does not converge.
- **xrd** -- Bragg's law, cubic d-spacings, and simulated powder-XRD
  patterns with systematic absences for SC/BCC/FCC, Laue interference,
  and Scherrer crystallite sizes.
- **defects** -- Frenkel and Schottky point-defect concentration vs.
  temperature.
- **crystal_chemistry** -- Pauling's radius-ratio rule and Goldschmidt's
  perovskite tolerance factor.
