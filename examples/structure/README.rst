Examples
========

This gallery walks through every public feature of
``chemistrykit.structure``: VSEPR geometry prediction with real 3D
coordinate generation; point-group determination from 3D coordinates
(genuine symmetry-element detection, not a formula lookup) and character
tables; bond order from the Pauling length correlation and from
Huckel-theory MO coefficients; formal-charge/oxidation-state
assignment from a Lewis structure and Pauling electronegativities;
Kekule-structure enumeration; Baeyer ring angle strain; and dipole
moments.

Each script in this gallery is self-contained and can be run directly with
``python examples/structure/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **vsepr** -- van't Hoff and Le Bel's tetrahedral carbon and chirality,
  steric number to idealized molecular shape, and hypervalent molecules.
- **point_group** -- Schoenflies symbols from detected symmetry elements,
  Werner's octahedral isomer counting, Mulliken's irrep labels, and
  Bethe's crystal-field splitting of the d orbitals.
- **bonding** -- the Coulson pi bond order from Huckel molecular-orbital
  theory, and the Pauling bond-order/bond-length correlation.
- **lewis** -- formal charges of competing Lewis structures, oxidation
  states in Kossel's ionic limit, and Pauling electronegativities from
  bond energies.
- **kekule** -- Kekule structures of benzene and larger aromatics.
- **ring_strain** -- Baeyer's angle strain in planar cycloalkanes.
- **dipole** -- molecular dipole moments from bond dipoles and charges.
