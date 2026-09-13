Examples
========

This gallery walks through every public feature of ``chemistrykit.thermo``:
ideal-gas, van der Waals, and Redlich-Kwong equations of state;
Clausius-Clapeyron phase boundaries; Gibbs-energy-minimization reaction
equilibrium; and Raoult's/Henry's law mixtures with colligative
properties.

Each script in this gallery is self-contained and can be run directly with
``python examples/thermo/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **equations_of_state** -- comparing the ideal gas law, van der Waals,
  and Redlich-Kwong P-V isotherms for a real gas, and the effect of
  approaching the critical point.
- **phase_equilibria** -- the Clausius-Clapeyron liquid-vapor phase
  boundary and the Gibbs phase rule.
- **equilibrium** -- the reaction quotient, Kp/Kc, the van't Hoff
  temperature dependence of an equilibrium constant, and the
  Gibbs-energy-minimization equilibrium-composition solver.
- **mixtures** -- Raoult's-law P-x-y diagrams for an ideal binary
  solution, the colligative properties (freezing-point depression,
  boiling-point elevation, osmotic pressure), and Henry's law vs.
  Raoult's law as the two limiting laws of a real solution.
