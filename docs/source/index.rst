chemistrykit
============

**chemistrykit** is a unified numerical toolkit for computational
chemistry, sharing common ODE integrators and chemical constants across
domain subpackages. This is an early, in-progress build -- see
``chemistrykit-spec.md`` in the repository root for the full 14-domain
plan; only the domains below exist so far.

- :mod:`chemistrykit.kinetics` -- reaction kinetics: integrated rate laws,
  the Arrhenius equation, Michaelis-Menten enzyme kinetics, a general
  stoichiometric reaction-network engine, and the Brusselator oscillator.
- :mod:`chemistrykit.thermo` -- chemical thermodynamics: equations of
  state (ideal gas, van der Waals, Redlich-Kwong), Clausius-Clapeyron
  phase boundaries, reaction equilibrium (Kp/Kc, van't Hoff, and a
  Gibbs-energy-minimization equilibrium-composition solver), and
  Raoult's/Henry's law mixtures with colligative properties.
- :mod:`chemistrykit.solutions` -- solution chemistry: pH/pOH and weak
  acid/base equilibria with Henderson-Hasselbalch buffers, titration
  curves, Ksp solubility equilibria, and Debye-Huckel activity
  coefficients.

Conventionally imported as ``ck``:

.. code-block:: python

   import chemistrykit as ck

   network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
   result = network.integrate((0.0, 10.0), dt=1e-3, method="rk4")

   acid = ck.solutions.WeakAcid(Ca=0.1, Ka=1.8e-5)
   print(acid.pH())

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/kinetics
   api/thermo
   api/solutions

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/kinetics
   examples/thermo
   examples/solutions
