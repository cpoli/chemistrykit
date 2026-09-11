"""chemistrykit.thermo: chemical thermodynamics.

Equations of state (ideal gas, van der Waals, Redlich-Kwong);
Clausius-Clapeyron phase boundaries and the Gibbs phase rule; reaction
equilibrium (Kp/Kc, the reaction quotient, the van't Hoff equation, and a
Gibbs-energy-minimization equilibrium-composition solver); and
Raoult's/Henry's law mixtures with colligative properties (freezing-point
depression, boiling-point elevation, osmotic pressure).
"""

__version__ = "0.1.0"

from chemistrykit.thermo.core.base_system import EquationOfState
from chemistrykit.thermo.systems.equations_of_state import IdealGas, RedlichKwong, VanDerWaals
from chemistrykit.thermo.systems.equilibrium import (
    EquilibriumComposition,
    VantHoffFit,
    fit_van_t_hoff,
    gibbs_energy_of_mixture,
    kc_from_kp,
    kp_from_kc,
    reaction_quotient,
    solve_equilibrium_composition,
    van_t_hoff_equilibrium_constant,
)
from chemistrykit.thermo.systems.mixtures import (
    CRYOSCOPIC_CONSTANTS,
    BinaryIdealSolution,
    boiling_point_elevation,
    freezing_point_depression,
    henry_law_pressure,
    osmotic_pressure,
    raoult_vapor_pressure,
)
from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron, gibbs_phase_rule

__all__ = [
    "__version__",
    "EquationOfState",
    "IdealGas",
    "VanDerWaals",
    "RedlichKwong",
    "ClausiusClapeyron",
    "gibbs_phase_rule",
    "reaction_quotient",
    "kp_from_kc",
    "kc_from_kp",
    "van_t_hoff_equilibrium_constant",
    "VantHoffFit",
    "fit_van_t_hoff",
    "EquilibriumComposition",
    "gibbs_energy_of_mixture",
    "solve_equilibrium_composition",
    "raoult_vapor_pressure",
    "BinaryIdealSolution",
    "henry_law_pressure",
    "CRYOSCOPIC_CONSTANTS",
    "freezing_point_depression",
    "boiling_point_elevation",
    "osmotic_pressure",
]
