"""chemistrykit: unified numerical toolkit for computational chemistry.

Import as ``ck`` by convention::

    import chemistrykit as ck
    ck.kinetics.FirstOrder(k=0.1, C0=1.0)
    ck.kinetics.StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
    ck.thermo.VanDerWaals(a=0.1448, b=3.913e-5)
    ck.solutions.WeakAcid(Ca=0.1, Ka=1.8e-5).pH()
    ck.md.LJFluid.from_lattice(n_per_side=4, density=0.6, temperature=1.0)
    ck.statmech.MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
    ck.structure.determine_point_group(water_molecule)
    ck.spectro.rotational_spectrum(rotor, J_max=10, temperature=300.0)
    ck.electrochem.nernst_potential(E_standard=0.34, n=2, Q=0.01)
    ck.photochem.jablonski_network(kf=2.0, kic=1.0, kisc=0.5, kp=0.3, kic_T=0.2)
    ck.surface.LangmuirIsotherm(K=2.0, qmax=5.0).loading(P=1.0)
    ck.polymer.IdealChain().end_to_end_distance(n=1000, b=0.5)
    ck.constants.R
    ck.integrators.rk4_integrate(...)

chemistrykit mirrors the architecture of the sibling project physicskit
(pk): one subpackage per chemistry domain, sharing common ODE integrators
(:mod:`chemistrykit.integrators`) and physical/chemical constants
(:mod:`chemistrykit.constants`). This is an early, in-progress build: only
the domains listed in ``__all__`` below exist so far -- see
``chemistrykit-spec.md`` for the full 14-domain plan and build order.
"""

from chemistrykit import (
    constants,
    electrochem,
    integrators,
    kinetics,
    md,
    photochem,
    polymer,
    quantum,
    solutions,
    spectro,
    statmech,
    structure,
    surface,
    thermo,
)

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
    "spectro",
    "structure",
    "electrochem",
    "photochem",
    "surface",
    "polymer",
]
