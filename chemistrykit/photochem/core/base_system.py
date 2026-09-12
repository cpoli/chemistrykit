"""Result dataclasses for chemistrykit.photochem; a note on why no new ABC is added here.

Every model in :mod:`chemistrykit.photochem` reduces, in the end, to
either a closed-form algebraic relationship (quantum yields, the
Stern-Volmer equation, the photostationary-state ratio) or a coupled
first-order rate-equation system of exactly the kind
:class:`chemistrykit.kinetics.core.base_system.ReactionNetwork` and
:class:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork`
already model generally -- a Jablonski-diagram excited-state decay is
literally a mass-action network with reactant order 1 throughout (see
:mod:`chemistrykit.photochem.systems.jablonski`), and a two-state
photoswitch under simultaneous forward/reverse photolysis is literally
:meth:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork.reversible`
with rate constants derived from photon-absorption/quantum-yield data
instead of Arrhenius kinetics (see
:mod:`chemistrykit.photochem.systems.photostationary_state`). Rather than
introduce a second, parallel network ABC purely for photochemistry, this
subpackage reuses ``chemistrykit.kinetics``'s directly -- following the
precedent of :mod:`chemistrykit.spectro` reusing
:mod:`chemistrykit.quantum`'s ``RigidRotor``/``QuantumHarmonicOscillator``
models rather than reimplementing quantum mechanics for spectroscopy.

The one genuinely new result type is
:class:`PhotostationaryStateResult`, a small dataclass (mirroring
:class:`chemistrykit.thermo.systems.equilibrium.EquilibriumComposition`)
bundling a photostationary-state calculation's populations, ratio, and
the rate constants that produced them.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["PhotostationaryStateResult"]


@dataclass
class PhotostationaryStateResult:
    """Result of a :func:`chemistrykit.photochem.systems.photostationary_state.photostationary_state` call."""

    A_pss: float
    """float: Steady-state concentration (or population fraction) of species A."""

    B_pss: float
    """float: Steady-state concentration (or population fraction) of species B."""

    ratio_B_over_A: float
    """float: :math:`[B]_{pss}/[A]_{pss} = k_{AB}/k_{BA}`."""

    k_AB: float
    """float: Forward (A -> B) pseudo-first-order photolysis rate constant used."""

    k_BA: float
    """float: Reverse (B -> A) pseudo-first-order photolysis rate constant used."""
