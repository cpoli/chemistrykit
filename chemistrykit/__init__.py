"""chemistrykit: unified numerical toolkit for computational chemistry.

Import as ``ck`` by convention::

    import chemistrykit as ck
    ck.constants.R
    ck.integrators.rk4_integrate(...)

chemistrykit mirrors the architecture of the sibling project physicskit
(pk): one subpackage per chemistry domain, sharing common ODE integrators
(:mod:`chemistrykit.integrators`) and physical/chemical constants
(:mod:`chemistrykit.constants`). This is an early, in-progress build: only
the domains listed in ``__all__`` below exist so far -- see
``chemistrykit-spec.md`` for the full 14-domain plan and build order.
"""

from chemistrykit import constants, integrators

__version__ = "0.1.0"

__all__ = [
    "constants",
    "integrators",
]
