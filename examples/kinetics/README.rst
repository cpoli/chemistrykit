Examples
========

This gallery walks through every public feature of ``chemistrykit.kinetics``:
closed-form integrated rate laws, the Arrhenius temperature dependence,
collision, diffusion, and transition-state theories of the rate constant,
Michaelis-Menten enzyme kinetics, a general stoichiometric reaction-network
engine with exact stochastic simulation, and chemical oscillators.

Each script in this gallery is self-contained and can be run directly with
``python examples/kinetics/<section>/<script>.py``. Every script also
carries an RST module docstring as its title/description and uses ``# %%``
markers to split narrative text from code, which is exactly what
Sphinx-Gallery renders into the pages below -- the script *is* the source of
truth for what you see, not a copy of it.

Sections
--------

- **rate_laws** -- the textbook zero/first/second-order integrated rate
  laws and their half-lives, reaction order from initial rates, and
  Wilhelmy's first-order sucrose inversion.
- **arrhenius** -- the Arrhenius temperature dependence of a rate constant,
  and recovering the activation energy from synthetic rate-vs-temperature
  data via an Arrhenius plot.
- **rate_theory** -- rate constants from molecular properties: collision
  theory, Smoluchowski's diffusion limit, and Eyring's transition-state
  theory.
- **enzyme** -- Michaelis-Menten enzyme kinetics, the Lineweaver-Burk
  linearization, competitive/noncompetitive inhibition, and the full
  substrate-depletion progress curve.
- **networks** -- the general stoichiometric reaction-network engine:
  the steady-state approximation and Lindemann fall-off, chain-branching
  explosions, Bateman's consecutive-reaction solution, Eigen's relaxation
  kinetics, parallel reactions, and Gillespie's stochastic simulation.
- **oscillators** -- the Brusselator limit cycle, Lotka's neutral
  oscillations, and the Oregonator model of the Belousov-Zhabotinsky
  reaction.
