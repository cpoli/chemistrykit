"""chemistrykit.polymer: polymer chemistry.

Ideal random-walk chain statistics (end-to-end distance, radius of
gyration) and the Flory exponent for real chains under different solvent
conditions; molecular-weight-distribution statistics (Mn, Mw, PDI) and
the Flory-Schulz (most-probable) distribution derived from step-growth
polymerization statistics; the Carothers equation for step-growth
kinetics; and chain-growth/free-radical polymerization kinetics built on
:mod:`chemistrykit.kinetics`'s reaction-network engine, with its
steady-state-approximation closed form as a cross-check.
"""

__version__ = "0.1.0"

from chemistrykit.polymer.core.base_system import PolymerChainModel
from chemistrykit.polymer.systems.chain_growth import (
    free_radical_network,
    kinetic_chain_length,
    steady_state_radical_concentration,
    steady_state_rate_of_polymerization,
)
from chemistrykit.polymer.systems.chain_statistics import FLORY_EXPONENTS, IdealChain, RealChain, flory_exponent
from chemistrykit.polymer.systems.molecular_weight_distribution import (
    MolecularWeightDistribution,
    flory_schulz_number_average_DP,
    flory_schulz_number_fraction,
    flory_schulz_pdi,
    flory_schulz_weight_average_DP,
    flory_schulz_weight_fraction,
    number_average_molar_mass,
    polydispersity_index,
    weight_average_molar_mass,
)
from chemistrykit.polymer.systems.step_growth import (
    degree_of_polymerization,
    degree_of_polymerization_stoichiometric_imbalance,
    extent_of_reaction_for_DP,
)

__all__ = [
    "__version__",
    "PolymerChainModel",
    "FLORY_EXPONENTS",
    "flory_exponent",
    "IdealChain",
    "RealChain",
    "number_average_molar_mass",
    "weight_average_molar_mass",
    "polydispersity_index",
    "MolecularWeightDistribution",
    "flory_schulz_number_fraction",
    "flory_schulz_weight_fraction",
    "flory_schulz_number_average_DP",
    "flory_schulz_weight_average_DP",
    "flory_schulz_pdi",
    "degree_of_polymerization",
    "extent_of_reaction_for_DP",
    "degree_of_polymerization_stoichiometric_imbalance",
    "free_radical_network",
    "steady_state_radical_concentration",
    "steady_state_rate_of_polymerization",
    "kinetic_chain_length",
]
