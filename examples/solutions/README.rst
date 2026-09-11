Examples
========

This gallery walks through every public feature of
``chemistrykit.solutions``: pH/pOH and weak acid/base equilibria with
Henderson-Hasselbalch buffers; strong/weak acid-base titration curves;
Ksp solubility equilibria and the common-ion effect; and Debye-Huckel
activity-coefficient laws.

Each script in this gallery is self-contained and can be run directly with
``python examples/solutions/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **acid_base** -- exact weak acid/base equilibrium (the Ostwald dilution
  law), and Henderson-Hasselbalch buffer design.
- **titration** -- strong/strong, weak-acid/strong-base, and
  weak-base/strong-acid titration curves, with numerically detected
  equivalence points.
- **solubility** -- Ksp, molar solubility, and the common-ion effect.
- **activity** -- Debye-Huckel limiting and extended activity-coefficient
  laws as a function of ionic strength.
