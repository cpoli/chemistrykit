"""chemistrykit.solutions: solution chemistry.

pH/pOH and weak acid/base equilibria with Henderson-Hasselbalch buffers,
Van Slyke buffer capacity, and Bjerrum polyprotic speciation; a full
titration-curve solver (strong/weak acid-base) with Gran-plot
equivalence-point estimation; solubility equilibria (Ksp, common-ion
effect); and Debye-Huckel limiting/extended and Davies
activity-coefficient laws. Colligative properties (freezing-point
depression, boiling-point elevation, osmotic pressure) live in
:mod:`chemistrykit.thermo` instead -- see
:func:`chemistrykit.thermo.freezing_point_depression` and friends.
"""

__version__ = "0.1.0"

from chemistrykit.solutions.core.base_system import GranPlotResult, Titration, TitrationResult, WeakElectrolyte
from chemistrykit.solutions.systems.acid_base import (
    Buffer,
    WeakAcid,
    WeakBase,
    buffer_capacity,
    h_from_ph,
    henderson_hasselbalch_ph,
    oh_from_poh,
    ph_from_h,
    poh_from_oh,
    polyprotic_fractions,
)
from chemistrykit.solutions.systems.activity import (
    DEBYE_HUCKEL_A_25C,
    activity_coefficient_davies,
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
    ionic_strength,
)
from chemistrykit.solutions.systems.solubility import (
    ksp_from_molar_solubility,
    molar_solubility_from_ksp,
    molar_solubility_with_common_ion,
)
from chemistrykit.solutions.systems.titration import (
    StrongAcidStrongBaseTitration,
    WeakAcidStrongBaseTitration,
    WeakBaseStrongAcidTitration,
    gran_plot,
)

__all__ = [
    "__version__",
    "WeakElectrolyte",
    "Titration",
    "TitrationResult",
    "GranPlotResult",
    "ph_from_h",
    "h_from_ph",
    "poh_from_oh",
    "oh_from_poh",
    "WeakAcid",
    "WeakBase",
    "Buffer",
    "henderson_hasselbalch_ph",
    "polyprotic_fractions",
    "buffer_capacity",
    "StrongAcidStrongBaseTitration",
    "WeakAcidStrongBaseTitration",
    "WeakBaseStrongAcidTitration",
    "gran_plot",
    "ksp_from_molar_solubility",
    "molar_solubility_from_ksp",
    "molar_solubility_with_common_ion",
    "DEBYE_HUCKEL_A_25C",
    "ionic_strength",
    "activity_coefficient_debye_huckel_limiting",
    "activity_coefficient_debye_huckel_extended",
    "activity_coefficient_davies",
]
