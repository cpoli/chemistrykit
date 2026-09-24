Examples
========

This gallery walks through every public feature of
``chemistrykit.analytical``: redox and complexometric (EDTA)
titration-curve simulation with equivalence-point detection and Gran
plots (acid-base titrations are covered by ``chemistrykit.solutions``'s
gallery); chromatographic plate theory, the van Deemter equation,
retention indices, and resolution; linear-regression calibration curves
with limits of detection/quantitation; propagation of uncertainty;
Dixon's and Grubbs' outlier tests; Student's t intervals and the Horwitz
function; and Savitzky-Golay smoothing.

Each script in this gallery is self-contained and can be run directly
with ``python examples/analytical/<section>/<script>.py``. Every script
also carries an RST module docstring as its title/description and uses
``# %%`` markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source
of truth for what you see, not a copy of it.

Sections
--------

- **titration** -- redox (Nernst) and EDTA titration curves with
  numerical equivalence-point detection, and the Gran plot.
- **chromatography** -- Tsvet's column separation, theoretical plates,
  the van Deemter and Golay equations, Kovats retention indices, and
  Purnell's resolution equation.
- **calibration** -- a least-squares calibration curve and its limits of
  detection/quantitation.
- **uncertainty** -- propagation-of-uncertainty formulas for sums,
  products, and powers, cross-checked against a general numerical
  implementation.
- **qtest** -- Dixon's Q-test and Grubbs' test for rejecting a suspect
  outlier from a small data set.
- **statistics** -- Student's t confidence intervals and the Horwitz
  interlaboratory-precision function.
- **smoothing** -- Savitzky-Golay smoothing and differentiation.
