Examples
========

This gallery walks through every public feature of ``chemistrykit.kinetics``:
closed-form integrated rate laws, the Arrhenius temperature dependence,
Michaelis-Menten enzyme kinetics, a general stoichiometric reaction-network
engine, and the Brusselator chemical oscillator.

Each script in this gallery is self-contained and can be run directly with
``python examples/kinetics/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **rate_laws** -- the textbook zero/first/second-order integrated rate
  laws and their half-lives, compared side by side.
- **arrhenius** -- the Arrhenius temperature dependence of a rate constant,
  and recovering the activation energy from synthetic rate-vs-temperature
  data via an Arrhenius plot.
- **enzyme** -- Michaelis-Menten enzyme kinetics, the Lineweaver-Burk
  linearization, competitive/noncompetitive inhibition, and the full
  substrate-depletion progress curve.
- **networks** -- the general stoichiometric reaction-network engine:
  parallel, consecutive (A->B->C), and reversible mechanisms, each checked
  against its closed-form solution, plus the steady-state approximation
  and a chain-branching explosion mechanism.
- **oscillators** -- the Brusselator: a chemical mechanism whose kinetics
  support a stable limit cycle rather than relaxing to equilibrium.
