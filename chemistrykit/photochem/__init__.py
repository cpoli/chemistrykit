"""chemistrykit.photochem: photochemistry.

Jablonski-diagram excited-state kinetics (fluorescence, internal
conversion, intersystem crossing, and phosphorescence) modeled as a
mass-action reaction network and integrated via
:mod:`chemistrykit.kinetics`/:mod:`chemistrykit.integrators`;
fluorescence/phosphorescence quantum yields and the photochemical
quantum yield via Beer-Lambert (:mod:`chemistrykit.spectro`);
Stern-Volmer fluorescence quenching with a static-vs-dynamic diagnostic;
photostationary-state kinetics for a two-state photoswitch under
simultaneous forward/reverse photolysis; the H2/Cl2 photochemical chain
reaction; Stokes shifts and Perrin fluorescence anisotropy; Kasha's-rule
emission yields; Förster and Dexter energy transfer; ferrioxalate
actinometry; and Rehm-Weller electron-transfer quenching.
"""

__version__ = "0.1.0"

from chemistrykit.photochem.core.base_system import PhotostationaryStateResult
from chemistrykit.photochem.systems.actinometry import ferrioxalate_fe2_moles, ferrioxalate_photon_flux
from chemistrykit.photochem.systems.chain_reaction import chain_quantum_yield, hydrogen_chlorine_chain_network
from chemistrykit.photochem.systems.electron_transfer import rehm_weller_free_energy, rehm_weller_quenching_rate
from chemistrykit.photochem.systems.energy_transfer import (
    dexter_rate,
    forster_efficiency,
    forster_radius,
    forster_rate,
)
from chemistrykit.photochem.systems.fluorescence import perrin_anisotropy, rotational_correlation_time, stokes_shift
from chemistrykit.photochem.systems.jablonski import (
    jablonski_network,
    jablonski_populations_analytic,
    kasha_emission_yields,
)
from chemistrykit.photochem.systems.photostationary_state import (
    photostationary_ratio,
    photostationary_state,
    photoswitch_network,
    photoswitch_rate_constants,
)
from chemistrykit.photochem.systems.quantum_yield import (
    fluorescence_quantum_yield,
    intersystem_crossing_yield,
    phosphorescence_quantum_yield,
    photochemical_quantum_yield,
    photons_absorbed,
)
from chemistrykit.photochem.systems.stern_volmer import (
    SternVolmerFit,
    classify_quenching_mechanism,
    dynamic_quenching_constant,
    fit_stern_volmer,
    stern_volmer_ratio,
)

__all__ = [
    "__version__",
    "PhotostationaryStateResult",
    "jablonski_network",
    "jablonski_populations_analytic",
    "kasha_emission_yields",
    "fluorescence_quantum_yield",
    "intersystem_crossing_yield",
    "phosphorescence_quantum_yield",
    "photons_absorbed",
    "photochemical_quantum_yield",
    "stern_volmer_ratio",
    "dynamic_quenching_constant",
    "SternVolmerFit",
    "fit_stern_volmer",
    "classify_quenching_mechanism",
    "photoswitch_rate_constants",
    "photoswitch_network",
    "photostationary_ratio",
    "photostationary_state",
    "hydrogen_chlorine_chain_network",
    "chain_quantum_yield",
    "stokes_shift",
    "perrin_anisotropy",
    "rotational_correlation_time",
    "forster_radius",
    "forster_rate",
    "forster_efficiency",
    "dexter_rate",
    "ferrioxalate_fe2_moles",
    "ferrioxalate_photon_flux",
    "rehm_weller_free_energy",
    "rehm_weller_quenching_rate",
]
