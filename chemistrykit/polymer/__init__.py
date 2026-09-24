"""chemistrykit.polymer: polymer chemistry.

Ideal random-walk chain statistics (end-to-end distance, radius of
gyration) and the Flory exponent for real chains under different solvent
conditions; molecular-weight-distribution statistics (Mn, Mw, PDI) and
the Flory-Schulz (most-probable) distribution derived from step-growth
polymerization statistics; the Carothers equation for step-growth
kinetics; and chain-growth/free-radical polymerization kinetics built on
:mod:`chemistrykit.kinetics`'s reaction-network engine, with its
steady-state-approximation closed form as a cross-check. Also Kuhn's
freely jointed chain and the Kratky-Porod worm-like chain; Flory-Huggins
solution thermodynamics; Staudinger/Mark-Houwink solution viscosity;
Debye light scattering vs. osmometry; Flory-Stockmayer gelation; the
Mayo-Lewis copolymer equation; Poisson (living) chain-length
distributions; and Bernoullian tacticity statistics.
"""

__version__ = "0.1.0"

from chemistrykit.polymer.core.base_system import PolymerChainModel
from chemistrykit.polymer.systems.chain_growth import (
    free_radical_network,
    kinetic_chain_length,
    steady_state_radical_concentration,
    steady_state_rate_of_polymerization,
)
from chemistrykit.polymer.systems.chain_statistics import (
    FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP,
    FLORY_EXPONENTS,
    IdealChain,
    RealChain,
    flory_exponent,
    freely_jointed_chain,
    worm_like_chain_mean_square_end_to_end,
)
from chemistrykit.polymer.systems.copolymerization import (
    azeotropic_feed_composition,
    mayo_lewis_copolymer_composition,
)
from chemistrykit.polymer.systems.flory_huggins import (
    flory_huggins_critical_point,
    flory_huggins_free_energy,
    flory_huggins_spinodal_chi,
)
from chemistrykit.polymer.systems.gelation import (
    branching_number_average_DP,
    branching_weight_average_DP,
    carothers_gel_point,
    flory_stockmayer_gel_point,
)
from chemistrykit.polymer.systems.light_scattering import (
    debye_Kc_over_R,
    osmotic_pressure_dilute_mixture,
    rayleigh_ratio_dilute_mixture,
)
from chemistrykit.polymer.systems.living_polymerization import (
    poisson_number_average_DP,
    poisson_number_fraction,
    poisson_pdi,
    poisson_weight_average_DP,
    simulate_living_polymerization,
)
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
from chemistrykit.polymer.systems.solution_viscosity import (
    fit_mark_houwink,
    mark_houwink_exponent_from_flory,
    mark_houwink_intrinsic_viscosity,
    staudinger_specific_viscosity,
)
from chemistrykit.polymer.systems.step_growth import (
    degree_of_polymerization,
    degree_of_polymerization_stoichiometric_imbalance,
    extent_of_reaction_for_DP,
)
from chemistrykit.polymer.systems.tacticity import (
    bernoullian_triad_fractions,
    mean_isotactic_run_length,
    sample_dyad_sequence,
)

__all__ = [
    "__version__",
    "PolymerChainModel",
    "FLORY_EXPONENTS",
    "FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP",
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
    "freely_jointed_chain",
    "worm_like_chain_mean_square_end_to_end",
    "flory_huggins_free_energy",
    "flory_huggins_spinodal_chi",
    "flory_huggins_critical_point",
    "staudinger_specific_viscosity",
    "mark_houwink_intrinsic_viscosity",
    "mark_houwink_exponent_from_flory",
    "fit_mark_houwink",
    "carothers_gel_point",
    "flory_stockmayer_gel_point",
    "branching_number_average_DP",
    "branching_weight_average_DP",
    "mayo_lewis_copolymer_composition",
    "azeotropic_feed_composition",
    "poisson_number_fraction",
    "poisson_number_average_DP",
    "poisson_weight_average_DP",
    "poisson_pdi",
    "simulate_living_polymerization",
    "bernoullian_triad_fractions",
    "mean_isotactic_run_length",
    "sample_dyad_sequence",
    "rayleigh_ratio_dilute_mixture",
    "debye_Kc_over_R",
    "osmotic_pressure_dilute_mixture",
]
