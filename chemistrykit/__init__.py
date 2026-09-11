"""chemistrykit: unified numerical toolkit for computational chemistry.

Import as ``ck`` by convention::

    import chemistrykit as ck
    ck.kinetics.FirstOrder(k=0.1, C0=1.0)
    ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
    ck.thermo.VanDerWaals(a=0.1448, b=3.913e-5)
    ck.solutions.WeakAcid(Ca=0.1, Ka=1.8e-5).pH()
    ck.md.LJFluid.from_lattice(n_per_side=4, density=0.6, temperature=1.0)
    ck.statmech.MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
    ck.constants.R
    ck.integrators.rk4_integrate(...)

chemistrykit mirrors the architecture of the sibling project physicskit
(pk): one subpackage per chemistry domain, sharing common ODE integrators
(:mod:`chemistrykit.integrators`) and physical/chemical constants
(:mod:`chemistrykit.constants`). This is an early, in-progress build: only
the domains listed in ``__all__`` below exist so far -- see
``chemistrykit-spec.md`` for the full 14-domain plan and build order.
"""

from chemistrykit import constants, integrators, kinetics, md, quantum, solutions, statmech, thermo

__version__ = "0.1.0"

__all__ = [
    "constants",
    "integrators",
    "kinetics",
    "thermo",
    "solutions",
    "md",
    "statmech",
    "quantum",
]
