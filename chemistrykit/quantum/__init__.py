"""chemistrykit.quantum: quantum chemistry -- pedagogical but numerically real.

Exactly solvable models (particle in a 1D/3D box, applied to the Kuhn
free-electron model of conjugated-dye UV-Vis absorption; the quantum
harmonic oscillator compared against the exact Morse-potential vibrational
levels; the rigid rotor; hydrogen-like radial wavefunctions and orbital
shapes) and genuinely variational ones (Huckel molecular-orbital theory
for conjugated pi systems, including Huckel's 4n+2 aromaticity rule
checked against the computed spectrum; a minimal Gaussian-basis
variational treatment of H2+; the effective-nuclear-charge variational
helium atom; Rayleigh-Schrodinger perturbation theory for
the anharmonic oscillator, checked against exact numerical diagonalization
in a truncated basis).

Every variational model here reduces to the same secular equation
:math:`HC=SCE` (:mod:`chemistrykit.quantum.utils.secular_equation`), which
:mod:`chemistrykit.quantum.systems.huckel`, :mod:`chemistrykit.quantum.systems.hartree_fock`,
and :mod:`chemistrykit.quantum.systems.perturbation` all build a matrix
Hamiltonian for and hand off to it -- mirroring the tight-binding
Hamiltonian-diagonalization pattern physicskit's
``physicskit.condensed.tight_binding`` uses for solid-state band
structure.
"""

__version__ = "0.1.0"

from chemistrykit.quantum.core.base_system import EigenstateResult, QuantumSystem, VariationalSolver
from chemistrykit.quantum.systems.harmonic_oscillator import MorseOscillator, QuantumHarmonicOscillator, compare_harmonic_vs_morse
from chemistrykit.quantum.systems.hartree_fock import ExponentOptimizationResult, H2PlusVariational
from chemistrykit.quantum.systems.helium import (
    HARTREE_ENERGY,
    HeliumVariationalResult,
    helium_like_variational_energy,
    optimize_helium_like_effective_charge,
)
from chemistrykit.quantum.systems.huckel import (
    HuckelSystem,
    cyclic_polyene_eigenvalues,
    is_aromatic_by_huckel_rule,
    linear_polyene_eigenvalues,
)
from chemistrykit.quantum.systems.hydrogenlike import HydrogenLikeAtom
from chemistrykit.quantum.systems.particle_in_box import ParticleInBox1D, ParticleInBox3D, conjugated_dye_absorption_wavelength
from chemistrykit.quantum.systems.perturbation import (
    anharmonic_energy_levels,
    anharmonic_hamiltonian_matrix,
    cubic_perturbation_first_order_correction,
    position_operator_matrix,
    quartic_perturbation_first_order_correction,
)
from chemistrykit.quantum.systems.rigid_rotor import RigidRotor
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation

__all__ = [
    "__version__",
    "QuantumSystem",
    "VariationalSolver",
    "EigenstateResult",
    "ParticleInBox1D",
    "ParticleInBox3D",
    "conjugated_dye_absorption_wavelength",
    "QuantumHarmonicOscillator",
    "MorseOscillator",
    "compare_harmonic_vs_morse",
    "RigidRotor",
    "HydrogenLikeAtom",
    "HuckelSystem",
    "linear_polyene_eigenvalues",
    "cyclic_polyene_eigenvalues",
    "is_aromatic_by_huckel_rule",
    "H2PlusVariational",
    "ExponentOptimizationResult",
    "HARTREE_ENERGY",
    "HeliumVariationalResult",
    "helium_like_variational_energy",
    "optimize_helium_like_effective_charge",
    "position_operator_matrix",
    "cubic_perturbation_first_order_correction",
    "quartic_perturbation_first_order_correction",
    "anharmonic_hamiltonian_matrix",
    "anharmonic_energy_levels",
    "solve_secular_equation",
]
