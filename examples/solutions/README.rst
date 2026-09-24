Examples
========

This gallery walks through every public feature of
``chemistrykit.solutions``: pH/pOH and weak acid/base equilibria with
Henderson-Hasselbalch buffers, buffer capacity, and polyprotic
speciation; strong/weak acid-base titration curves and Gran plots; Ksp
solubility equilibria and the common-ion effect; and Debye-Huckel and
Davies activity-coefficient laws.

Each script in this gallery is self-contained and can be run directly with
``python examples/solutions/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **acid_base** -- mass action, Arrhenius dissociation, the Ostwald
  dilution law, the pH scale, conjugate pairs, Henderson-Hasselbalch
  buffers, buffer capacity, and polyprotic speciation.
- **titration** -- titration curves with numerically detected
  equivalence points, and Gran-plot extrapolation.
- **solubility** -- the solubility product and the common-ion effect.
- **activity** -- Debye-Huckel limiting, Guntelberg extended, and Davies
  activity-coefficient laws as a function of ionic strength.
