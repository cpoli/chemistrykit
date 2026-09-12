Examples
========

This gallery walks through every public feature of
``chemistrykit.structure``: VSEPR geometry prediction with real 3D
coordinate generation; point-group determination from 3D coordinates
(genuine symmetry-element detection, not a formula lookup) and character
tables; bond order from the Pauling length correlation and from
Huckel-theory MO coefficients; and formal-charge/oxidation-state
assignment from a Lewis structure.

Each script in this gallery is self-contained and can be run directly with
``python examples/structure/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **vsepr** -- steric number to idealized molecular shape, with real 3D
  coordinates for methane, water, sulfur tetrafluoride, and xenon
  tetrafluoride.
- **point_group** -- symmetry-element detection and point-group
  classification for water, ammonia, methane, and carbon dioxide.
- **bonding** -- the Pauling bond-order/bond-length correlation, and the
  Coulson pi bond order from Huckel molecular-orbital theory.
- **lewis** -- formal charge and oxidation-state assignment for a handful
  of classic Lewis-structure examples.
