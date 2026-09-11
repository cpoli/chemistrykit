chemistrykit
============

**chemistrykit** is a unified numerical toolkit for computational
chemistry, sharing common ODE integrators and chemical constants across
domain subpackages. This is an early, in-progress build -- see
``chemistrykit-spec.md`` in the repository root for the full 14-domain
plan; only the domain below exists so far.

- :mod:`chemistrykit.kinetics` -- reaction kinetics: integrated rate laws,
  the Arrhenius equation, Michaelis-Menten enzyme kinetics, a general
  stoichiometric reaction-network engine, and the Brusselator oscillator.

Conventionally imported as ``ck``:

.. code-block:: python

   import chemistrykit as ck

   network = ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
   result = network.integrate((0.0, 10.0), dt=1e-3, method="rk4")

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/kinetics

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/kinetics
