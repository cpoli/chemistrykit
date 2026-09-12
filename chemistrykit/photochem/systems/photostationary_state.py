r"""Photostationary-state kinetics for a two-state photoswitch under simultaneous forward/reverse photolysis.

A photoswitch (e.g. an azobenzene or diarylethene) interconverts between
two forms, :math:`A` and :math:`B`, each of which can absorb the
excitation light and photoisomerize to the other -- unlike a simple
photoreaction that runs to completion, both directions run
*simultaneously* under continuous illumination, so the system reaches a
photostationary state (PSS), a dynamic (not thermodynamic) equilibrium
with a nonzero net photon flux driving nonzero forward and reverse
photoisomerization rates that happen to exactly cancel in the
population balance (Turro, Ramamurthy & Scaiano, *Modern Molecular
Photochemistry of Organic Molecules* (2010), Ch. 7.6; E. Fischer, *J.
Phys. Chem.* 71, 3704 (1967), for the composition formula used here).

In the low-optical-density (dilute) limit -- every absorbed photon
converted with quantum yield :math:`\Phi`, and self-absorption between
the two species neglected -- each direction is a pseudo-first-order
process with rate constant :math:`k=2.303\,\Phi\,\varepsilon\,I_0` (`I_0`
the incident photon flux, `\varepsilon` the molar absorptivity at the
irradiation wavelength), so :math:`A\rightleftharpoons B` is *exactly*
the same math as
:meth:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork.reversible`
-- this module builds that network directly and derives the PSS ratio
from its known equilibrium-population formula rather than reimplementing
reversible first-order kinetics.

**Approximation flagged**: the pseudo-first-order rate constants below
are only strictly valid in the low-optical-density limit; at higher
concentration the two species compete for the same absorbed photons in a
way that makes the "rate constants" themselves concentration-dependent
(a more elaborate treatment than implemented here -- see Turro et al.,
Ch. 7.6, for the general inner-filter-corrected case).
"""

from __future__ import annotations

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork, reversible_analytic
from chemistrykit.photochem.core.base_system import PhotostationaryStateResult

__all__ = ["photoswitch_rate_constants", "photoswitch_network", "photostationary_ratio", "photostationary_state"]


def photoswitch_rate_constants(phi_AB: float, eps_A: float, phi_BA: float, eps_B: float, I0: float = 1.0) -> tuple[float, float]:
    r"""Pseudo-first-order photoisomerization rate constants in the low-optical-density limit.

    :math:`k_{AB}=2.303\,\Phi_{AB}\,\varepsilon_A\,I_0` and
    :math:`k_{BA}=2.303\,\Phi_{BA}\,\varepsilon_B\,I_0` (Turro, Ramamurthy
    & Scaiano, *Modern Molecular Photochemistry of Organic Molecules*,
    Ch. 7.6) -- see the module docstring for the dilute-limit
    approximation this relies on.

    Parameters
    ----------
    phi_AB : float
        Quantum yield of the :math:`A\to B` photoreaction.
    eps_A : float
        Molar absorptivity of :math:`A` at the irradiation wavelength,
        in L mol^-1 cm^-1.
    phi_BA : float
        Quantum yield of the :math:`B\to A` photoreaction.
    eps_B : float
        Molar absorptivity of :math:`B` at the irradiation wavelength.
    I0 : float, default 1.0
        Incident photon flux (or intensity); the ratio ``k_AB/k_BA`` --
        and hence the photostationary ratio -- does not depend on this
        value.

    Returns
    -------
    k_AB, k_BA : float

    Examples
    --------
    Equal quantum yields and absorptivities give equal rate constants:

    >>> k_AB, k_BA = photoswitch_rate_constants(phi_AB=0.5, eps_A=1000.0, phi_BA=0.5, eps_B=1000.0)
    >>> round(k_AB, 6) == round(k_BA, 6)
    True
    """
    k_AB = 2.303 * phi_AB * eps_A * I0
    k_BA = 2.303 * phi_BA * eps_B * I0
    return k_AB, k_BA


def photoswitch_network(k_AB: float, k_BA: float, A0: float = 1.0) -> StoichiometricNetwork:
    r"""Build the :math:`A\rightleftharpoons B` photoswitch network as a :class:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork`.

    Directly reuses
    :meth:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork.reversible`
    -- a two-state photoswitch under simultaneous forward/reverse
    photolysis is, mathematically, the identical reversible first-order
    reaction network that thermal isomerization would be, just with
    photochemical rather than Arrhenius rate constants.

    Parameters
    ----------
    k_AB, k_BA : float
        Pseudo-first-order photoisomerization rate constants (see
        :func:`photoswitch_rate_constants`).
    A0 : float, default 1.0
        Total (initial, all-A) population/concentration.

    Returns
    -------
    StoichiometricNetwork

    Examples
    --------
    >>> net = photoswitch_network(k_AB=1.5, k_BA=0.4, A0=1.0)
    >>> result = net.integrate((0.0, 50.0), dt=1e-3, method="rk4")
    >>> round(float(result.concentration("B")[-1] / result.concentration("A")[-1]), 3)
    3.75
    """
    return StoichiometricNetwork.reversible(kf=k_AB, kr=k_BA, A0=A0, B0=0.0)


def photostationary_ratio(k_AB: float, k_BA: float) -> float:
    r"""The photostationary-state population ratio :math:`[B]_{pss}/[A]_{pss} = k_{AB}/k_{BA}`.

    At the photostationary state, :math:`d[A]/dt=-k_{AB}[A]+k_{BA}[B]=0`,
    giving this ratio directly (E. Fischer, *J. Phys. Chem.* 71, 3704
    (1967); Turro, Ramamurthy & Scaiano, *Modern Molecular Photochemistry
    of Organic Molecules*, Ch. 7.6) -- exactly the long-time limit of
    :func:`chemistrykit.kinetics.systems.networks.reversible_analytic`'s
    equilibrium-population formula, since a photostationary state is
    mathematically the same steady state a thermally reversible
    first-order reaction relaxes to.

    Parameters
    ----------
    k_AB, k_BA : float
        Pseudo-first-order photoisomerization rate constants.

    Returns
    -------
    float

    Examples
    --------
    Equal rate constants give a 1:1 photostationary mixture:

    >>> round(photostationary_ratio(k_AB=2.0, k_BA=2.0), 6)
    1.0
    """
    return k_AB / k_BA


def photostationary_state(k_AB: float, k_BA: float, total_concentration: float = 1.0) -> PhotostationaryStateResult:
    r"""Solve for the photostationary-state populations of a two-state photoswitch.

    Uses the exact algebraic steady-state solution of
    :math:`A\rightleftharpoons B` (the :math:`t\to\infty` limit of
    :func:`chemistrykit.kinetics.systems.networks.reversible_analytic`,
    with the "forward" rate constant :math:`k_f=k_{AB}` and "reverse"
    :math:`k_r=k_{BA}`) rather than integrating :func:`photoswitch_network`
    to long time numerically -- both were checked to agree (see this
    module's test suite) as part of developing this function.

    Parameters
    ----------
    k_AB, k_BA : float
        Pseudo-first-order photoisomerization rate constants (see
        :func:`photoswitch_rate_constants`).
    total_concentration : float, default 1.0
        Total (conserved) :math:`[A]+[B]`.

    Returns
    -------
    PhotostationaryStateResult

    Examples
    --------
    >>> result = photostationary_state(k_AB=1.5, k_BA=0.4, total_concentration=1.0)
    >>> round(result.ratio_B_over_A, 3)
    3.75
    >>> round(result.A_pss + result.B_pss, 9)
    1.0
    """
    A_pss, B_pss = reversible_analytic(A0=total_concentration, kf=k_AB, kr=k_BA, t=1.0e18, B0=0.0)
    return PhotostationaryStateResult(
        A_pss=float(A_pss),
        B_pss=float(B_pss),
        ratio_B_over_A=photostationary_ratio(k_AB, k_BA),
        k_AB=k_AB,
        k_BA=k_BA,
    )
