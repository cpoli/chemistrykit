chemistrykit
============

**chemistrykit** is a unified numerical toolkit for computational
chemistry, sharing common ODE integrators and chemical constants across
domain subpackages. This is an early, in-progress build -- see
``chemistrykit-spec.md`` in the repository root for the full 14-domain
plan; no domain subpackage has landed yet.

Conventionally imported as ``ck``:

.. code-block:: python

   import chemistrykit as ck

   ck.constants.R
   ck.integrators.rk4_integrate(...)
