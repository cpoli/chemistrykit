"""chemistrykit.statmech: statistical mechanics of molecules.

Translational/rotational/vibrational partition functions and the
thermodynamic functions (U, S, Cv) derived from them (including the
Sackur-Tetrode translational entropy); the Maxwell-Boltzmann speed
distribution; and a canonical-ensemble lattice-gas adsorption model
(a statistical-mechanical derivation of the Langmuir isotherm).

:mod:`chemistrykit.md`'s trajectories are a natural cross-check for this
domain's distribution functions -- see
``examples/md/lj_fluid/plot_02_maxwell_boltzmann_check.py``.
"""

__version__ = "0.1.0"

from chemistrykit.statmech.core.base_system import PartitionFunction, ThermodynamicFunctions
from chemistrykit.statmech.systems.lattice_gas import LatticeGasAdsorption
from chemistrykit.statmech.systems.maxwell_boltzmann import (
    MaxwellBoltzmannSpeedDistribution,
    mean_speed,
    most_probable_speed,
    rms_speed,
)
from chemistrykit.statmech.systems.partition_functions import (
    IdealGasMolecule,
    RotationalPartitionFunctionLinear,
    TranslationalPartitionFunction,
    VibrationalPartitionFunctionHarmonic,
    sackur_tetrode_entropy,
)
from chemistrykit.statmech.utils.combinatorics import ln_binomial, ln_factorial
from chemistrykit.statmech.utils.thermal_wavelength import thermal_de_broglie_wavelength

__all__ = [
    "__version__",
    "PartitionFunction",
    "ThermodynamicFunctions",
    "TranslationalPartitionFunction",
    "RotationalPartitionFunctionLinear",
    "VibrationalPartitionFunctionHarmonic",
    "IdealGasMolecule",
    "sackur_tetrode_entropy",
    "MaxwellBoltzmannSpeedDistribution",
    "most_probable_speed",
    "mean_speed",
    "rms_speed",
    "LatticeGasAdsorption",
    "thermal_de_broglie_wavelength",
    "ln_factorial",
    "ln_binomial",
]
