"""chemistrykit.electrochem: electrochemistry.

The Nernst equation for standard and concentration cells (with an
activity-coefficient-corrected variant built on
:mod:`chemistrykit.solutions.systems.activity`); a curated table of
standard reduction potentials with redox-couple balancing and
cell-potential combination; Butler-Volmer electrode kinetics, exchange
current density, and Tafel-plot linearization; galvanic vs. electrolytic
cells and Faraday's laws of electrolysis; and a simplified constant-
current battery discharge/capacity model with Peukert's-law rate
dependence; Kohlrausch's conductivity laws; fuel-cell thermodynamic
limits; and the diffusion-limited currents of electroanalysis (Cottrell,
Ilkovič, Heyrovský-Ilkovič wave, Randles-Ševčík).
"""

__version__ = "0.1.0"

from chemistrykit.electrochem.core.base_system import BatteryDischargeModel, DischargeResult
from chemistrykit.electrochem.systems.battery import ConstantCurrentBattery, effective_capacity, peukert_discharge_time
from chemistrykit.electrochem.systems.butler_volmer import (
    TafelFit,
    butler_volmer_current_density,
    exchange_current_density,
    fit_tafel_plot,
    tafel_overpotential,
    tafel_slope,
)
from chemistrykit.electrochem.systems.conductivity import (
    LIMITING_IONIC_CONDUCTIVITIES,
    KohlrauschFit,
    fit_kohlrausch_law,
    kohlrausch_molar_conductivity,
    limiting_molar_conductivity,
)
from chemistrykit.electrochem.systems.electrolysis import (
    charge_from_current,
    faradays_law_mass,
    mass_from_charge,
    minimum_applied_voltage_electrolytic,
    moles_from_charge,
)
from chemistrykit.electrochem.systems.fuel_cell import (
    fuel_cell_efficiency_limit,
    reversible_cell_voltage,
    reversible_cell_voltage_at_temperature,
)
from chemistrykit.electrochem.systems.nernst import (
    activity_corrected_reaction_quotient,
    concentration_cell_potential,
    nernst_potential,
    nernst_potential_with_activity,
)
from chemistrykit.electrochem.systems.standard_potentials import (
    STANDARD_REDUCTION_POTENTIALS,
    HalfReaction,
    balance_redox_reaction,
    cell_potential,
    is_spontaneous,
    standard_cell_potential,
)
from chemistrykit.electrochem.systems.voltammetry import (
    cottrell_current,
    ilkovic_diffusion_current,
    polarographic_wave_current,
    randles_sevcik_peak_current,
)

__all__ = [
    "__version__",
    "BatteryDischargeModel",
    "DischargeResult",
    "nernst_potential",
    "concentration_cell_potential",
    "activity_corrected_reaction_quotient",
    "nernst_potential_with_activity",
    "HalfReaction",
    "STANDARD_REDUCTION_POTENTIALS",
    "cell_potential",
    "standard_cell_potential",
    "balance_redox_reaction",
    "is_spontaneous",
    "butler_volmer_current_density",
    "exchange_current_density",
    "tafel_slope",
    "tafel_overpotential",
    "TafelFit",
    "fit_tafel_plot",
    "charge_from_current",
    "moles_from_charge",
    "mass_from_charge",
    "faradays_law_mass",
    "minimum_applied_voltage_electrolytic",
    "peukert_discharge_time",
    "effective_capacity",
    "ConstantCurrentBattery",
    "LIMITING_IONIC_CONDUCTIVITIES",
    "limiting_molar_conductivity",
    "kohlrausch_molar_conductivity",
    "KohlrauschFit",
    "fit_kohlrausch_law",
    "reversible_cell_voltage",
    "reversible_cell_voltage_at_temperature",
    "fuel_cell_efficiency_limit",
    "cottrell_current",
    "ilkovic_diffusion_current",
    "polarographic_wave_current",
    "randles_sevcik_peak_current",
]
