"""chemistrykit.kinetics: reaction kinetics.

Elementary integrated rate laws and half-lives, the Arrhenius equation and
activation-energy fitting, Michaelis-Menten enzyme kinetics with
Lineweaver-Burk linearization and inhibition, a general
stoichiometry-matrix reaction-network engine (parallel, consecutive,
reversible, and steady-state-approximation chain reactions) integrated
via :mod:`chemistrykit.integrators`, and the Brusselator oscillating
reaction network.
"""

__version__ = "0.1.0"

from chemistrykit.kinetics.core.base_system import KineticsResult, RateLaw, ReactionNetwork
from chemistrykit.kinetics.systems.arrhenius import ArrheniusFit, arrhenius_rate_constant, fit_arrhenius
from chemistrykit.kinetics.systems.enzyme import (
    MichaelisMentenFit,
    MichaelisMentenProgress,
    competitive_inhibition_rate,
    fit_lineweaver_burk,
    michaelis_menten_rate,
    noncompetitive_inhibition_rate,
)
from chemistrykit.kinetics.systems.networks import (
    StoichiometricNetwork,
    consecutive_analytic,
    reversible_analytic,
    ssa_intermediate_concentration,
)
from chemistrykit.kinetics.systems.oscillators import Brusselator
from chemistrykit.kinetics.systems.rate_laws import FirstOrder, SecondOrder, ZeroOrder

__all__ = [
    "__version__",
    "KineticsResult",
    "RateLaw",
    "ReactionNetwork",
    "ZeroOrder",
    "FirstOrder",
    "SecondOrder",
    "arrhenius_rate_constant",
    "ArrheniusFit",
    "fit_arrhenius",
    "michaelis_menten_rate",
    "competitive_inhibition_rate",
    "noncompetitive_inhibition_rate",
    "MichaelisMentenFit",
    "fit_lineweaver_burk",
    "MichaelisMentenProgress",
    "StoichiometricNetwork",
    "consecutive_analytic",
    "reversible_analytic",
    "ssa_intermediate_concentration",
    "Brusselator",
]
