Examples
========

This gallery walks through every public feature of
``chemistrykit.analytical``: redox and complexometric (EDTA)
titration-curve simulation with equivalence-point detection (acid-base
titrations are covered by ``chemistrykit.solutions``'s gallery, and used
here directly alongside the new titration types); chromatographic plate
theory and the van Deemter equation; linear-regression calibration
curves with IUPAC-convention limits of detection/quantitation; and
propagation of uncertainty plus Dixon's Q-test for outlier rejection.

Each script in this gallery is self-contained and can be run directly
with ``python examples/analytical/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **titration** -- acid-base, redox, and EDTA titration curves, side by
  side, with numerical equivalence-point detection.
- **chromatography** -- theoretical plates, the van Deemter equation and
  its optimum flow velocity, and resolution/selectivity between two peaks.
- **calibration** -- a least-squares calibration curve and its IUPAC
  limits of detection/quantitation.
- **uncertainty** -- propagation-of-uncertainty formulas for sums,
  products, and powers, cross-checked against a general numerical
  implementation.
- **qtest** -- Dixon's Q-test for rejecting a suspect outlier from a
  small data set.
