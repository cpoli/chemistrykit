"""chemistrykit.md: molecular dynamics and force fields.

Lennard-Jones fluid simulation in reduced units (energy, pressure, radial
distribution function g(r)); pairwise potentials (Morse, Buckingham,
harmonic bond/angle terms); SHAKE bond constraints; velocity-rescaling,
Berendsen, Nose-Hoover, and stochastic velocity-rescaling thermostats;
Einstein and Green-Kubo self-diffusion coefficients; periodic boundary conditions with the minimum-image
convention and a Verlet neighbor list. Integrated via
:mod:`chemistrykit.integrators`'s velocity-Verlet
(:func:`chemistrykit.integrators.velocity_verlet_step`) -- not a
reimplementation, see
:meth:`chemistrykit.md.core.base_system.MolecularDynamicsSystem.step`.

MD's trajectories are a natural cross-check for
:mod:`chemistrykit.statmech`'s distribution functions: the speed
distribution of an equilibrated :class:`~chemistrykit.md.systems.lj_fluid.LJFluid`
should follow :class:`chemistrykit.statmech.MaxwellBoltzmannSpeedDistribution`
at the simulation's instantaneous temperature (see
``examples/md/lj_fluid/plot_02_maxwell_boltzmann_check.py``).
"""

__version__ = "0.1.0"

from chemistrykit.md.core.base_system import AnglePotential, MDResult, MolecularDynamicsSystem, PairPotential
from chemistrykit.md.systems.constraints import ShakeMolecule, shake
from chemistrykit.md.systems.lj_fluid import LennardJones, LJFluid
from chemistrykit.md.systems.pair_potentials import (
    Buckingham,
    DiatomicOscillator,
    HarmonicAngle,
    HarmonicBond,
    HarmonicMolecule,
    Morse,
)
from chemistrykit.md.systems.thermostats import (
    BerendsenThermostat,
    NoseHooverThermostat,
    StochasticVelocityRescalingThermostat,
    VelocityRescalingThermostat,
)
from chemistrykit.md.systems.transport import (
    einstein_diffusion_coefficient,
    green_kubo_diffusion_coefficient,
    mean_squared_displacement,
    unwrap_trajectory,
    velocity_autocorrelation,
)
from chemistrykit.md.utils.neighbor_list import VerletNeighborList, build_neighbor_list
from chemistrykit.md.utils.pbc import minimum_image_displacement, wrap_positions

__all__ = [
    "__version__",
    "MDResult",
    "PairPotential",
    "AnglePotential",
    "MolecularDynamicsSystem",
    "LennardJones",
    "LJFluid",
    "Morse",
    "Buckingham",
    "HarmonicBond",
    "HarmonicAngle",
    "DiatomicOscillator",
    "HarmonicMolecule",
    "VelocityRescalingThermostat",
    "NoseHooverThermostat",
    "BerendsenThermostat",
    "StochasticVelocityRescalingThermostat",
    "shake",
    "ShakeMolecule",
    "unwrap_trajectory",
    "mean_squared_displacement",
    "velocity_autocorrelation",
    "einstein_diffusion_coefficient",
    "green_kubo_diffusion_coefficient",
    "build_neighbor_list",
    "VerletNeighborList",
    "minimum_image_displacement",
    "wrap_positions",
]
