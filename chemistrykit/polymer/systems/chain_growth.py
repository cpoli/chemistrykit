r"""Chain-growth (free-radical) polymerization kinetics, and its steady-state approximation.

See Odian, *Principles of Polymerization*, 4th ed., Ch. 3.2-3.3
("Kinetics of chain polymerization"). Free-radical polymerization
proceeds through three elementary steps:

* **Initiation**: an initiator :math:`I` decomposes (rate constant
  :math:`k_d`) into two primary radicals, of which a fraction :math:`f`
  (the initiator efficiency; the rest are lost to side reactions such as
  cage recombination) go on to add a monomer and start a growing chain.
  The overall rate of radical *production* is
  :math:`R_i=2fk_d[I]`.
* **Propagation**: a growing radical chain :math:`M_n^\bullet` adds
  another monomer, :math:`M_n^\bullet+M\to M_{n+1}^\bullet`, at rate
  :math:`k_p[M^\bullet][M]` -- with :math:`k_p` taken independent of
  chain length `n` (the standard "long-chain" approximation), every
  propagation step looks kinetically identical, so the *total* radical
  concentration :math:`[M^\bullet]=\sum_n[M_n^\bullet]` is the only
  radical-related quantity that matters for the monomer/radical/polymer
  concentration evolution modeled here (individual chain lengths are not
  tracked -- see :mod:`chemistrykit.polymer.systems.molecular_weight_distribution`
  for the resulting distribution's closed form, in the step-growth case).
* **Termination**: two radical chains combine or disproportionate into
  dead polymer at rate :math:`k_t[M^\bullet]^2`, consuming radicals at
  rate :math:`2k_t[M^\bullet]^2`.

:func:`free_radical_network` builds this three-species-plus-dead-polymer
system as a
:class:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork` and
integrates it numerically, reusing the kinetics domain's general
mass-action reaction-network engine rather than hand-rolling new ODE
wrapping (per this package's standing convention -- see
:mod:`chemistrykit.photochem.systems.jablonski` for the precedent).
Because radical decomposition (:math:`k_d`) is typically many orders of
magnitude slower than termination (:math:`k_t[M^\bullet]`), the radical
concentration relaxes to a slowly-drifting quasi-steady-state almost
immediately (:math:`d[M^\bullet]/dt\approx0`, the **steady-state
approximation**, SSA): setting :math:`R_i=R_t` gives the closed-form
:func:`steady_state_radical_concentration`, used in
:mod:`chemistrykit.polymer.tests.test_chain_growth` as the cross-check
against direct numerical integration of :func:`free_radical_network`
-- exactly the same style of SSA-vs-numerical-integration check
:func:`chemistrykit.kinetics.systems.networks.ssa_intermediate_concentration`
performs for a simple consecutive reaction chain.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork

__all__ = [
    "free_radical_network",
    "steady_state_radical_concentration",
    "steady_state_rate_of_polymerization",
    "kinetic_chain_length",
]


def free_radical_network(kd: float, f: float, kp: float, kt: float, I0: float, M0: float, mode: str = "combination") -> StoichiometricNetwork:
    r"""Build the lumped initiation/propagation/termination free-radical polymerization network.

    Species, in order: ``I`` (initiator), ``M`` (monomer), ``R`` (total
    radical concentration, lumped over all chain lengths), ``D`` (dead
    polymer chains formed by termination -- a count of *chains*, not
    monomer units incorporated). The three reactions:

    1. ``I -> 2R`` at rate :math:`(fk_d)[I]` (so the net radical
       production rate is :math:`2fk_d[I]=R_i`, matching the standard
       definition above).
    2. ``R + M -> R`` at rate :math:`k_p[R][M]` (propagation: consumes a
       monomer, leaves the radical count unchanged -- the radical that
       reacted is regenerated one unit longer).
    3. ``2R -> D`` (`mode="combination"`, one dead chain per termination
       event) or ``2R -> 2D`` (`mode="disproportionation"`, two dead
       chains per event) at rate :math:`k_t[R]^2`.

    Parameters
    ----------
    kd : float
        Initiator decomposition rate constant.
    f : float
        Initiator efficiency, :math:`0<f\le1`.
    kp : float
        Propagation rate constant.
    kt : float
        Termination rate constant.
    I0, M0 : float
        Initial initiator and monomer concentrations.
    mode : {"combination", "disproportionation"}
        Termination mechanism (affects only the ``D`` stoichiometry, not
        the radical/monomer dynamics).

    Returns
    -------
    StoichiometricNetwork
        Species ``("I", "M", "R", "D")``.

    Examples
    --------
    Radical concentration quickly rises from zero and plateaus near the
    steady-state approximation's prediction (see
    :func:`steady_state_radical_concentration`), while initiator and
    monomer are consumed only slowly on that timescale:

    >>> kd, f, kp, kt = 1.0e-5, 0.5, 1.0e3, 1.0e7
    >>> I0, M0 = 0.01, 5.0
    >>> net = free_radical_network(kd, f, kp, kt, I0, M0)
    >>> result = net.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    >>> R_ss = steady_state_radical_concentration(kd, f, I0, kt)
    >>> bool(np.isclose(result.concentration("R")[-1], R_ss, rtol=0.05))
    True
    >>> bool(result.concentration("I")[-1] > 0.99 * I0)
    True
    """
    if mode not in ("combination", "disproportionation"):
        raise ValueError("mode must be 'combination' or 'disproportionation'")
    species = ("I", "M", "R", "D")
    D_per_event = 1.0 if mode == "combination" else 2.0
    stoich = [
        [-1.0, 0.0, 0.0],  # I
        [0.0, -1.0, 0.0],  # M
        [2.0, 0.0, -2.0],  # R
        [0.0, 0.0, D_per_event],  # D
    ]
    reactant_orders = [
        [1.0, 0.0, 0.0],  # I: order 1 in reaction 1
        [0.0, 1.0, 0.0],  # M: order 1 in reaction 2
        [0.0, 1.0, 2.0],  # R: order 1 in reaction 2, order 2 in reaction 3
        [0.0, 0.0, 0.0],  # D: never a reactant
    ]
    rate_constants = [f * kd, kp, kt]
    state0 = [I0, M0, 0.0, 0.0]
    return StoichiometricNetwork(species, stoich, rate_constants, reactant_orders, state0)


def steady_state_radical_concentration(kd: float, f: float, I, kt: float):
    r"""Steady-state-approximation (SSA) radical concentration :math:`[M^\bullet]_{ss}=\sqrt{fk_d[I]/k_t}`.

    Setting the radical production rate :math:`R_i=2fk_d[I]` equal to the
    termination (radical-consumption) rate :math:`R_t=2k_t[M^\bullet]^2`
    and solving for :math:`[M^\bullet]` gives this result (Odian,
    *Principles of Polymerization*, 4th ed., Ch. 3.3.1, eq. 3-21) -- valid
    whenever radical production and consumption equilibrate on a
    timescale much shorter than the initiator/monomer are consumed
    (:math:`k_t[M^\bullet]\gg k_d`, essentially always true in practice
    since termination is diffusion-controlled while initiator
    decomposition is not).

    Parameters
    ----------
    kd : float
        Initiator decomposition rate constant.
    f : float
        Initiator efficiency.
    I : float or array-like of float
        Initiator concentration.
    kt : float
        Termination rate constant.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(steady_state_radical_concentration(kd=1.0e-5, f=0.5, I=0.01, kt=1.0e7)), 10)
    7.07e-08
    """
    I = np.asarray(I, dtype=np.float64)
    result = np.sqrt(f * kd * I / kt)
    return float(result) if result.ndim == 0 else result


def steady_state_rate_of_polymerization(kd: float, f: float, kp: float, kt: float, I, M):
    r"""Steady-state rate of polymerization :math:`R_p=k_p[M]\sqrt{fk_d[I]/k_t}`.

    The rate at which monomer is consumed by propagation, evaluated at
    the steady-state radical concentration
    (:func:`steady_state_radical_concentration`) -- the classic result
    that :math:`R_p` scales as the *square root* of initiator
    concentration (Odian, *Principles of Polymerization*, 4th ed., Ch.
    3.3.1, eq. 3-22), a distinctive experimental signature of free-radical
    (as opposed to, e.g., ionic) chain polymerization.

    Parameters
    ----------
    kd, f, kp, kt : float
        Rate constants and initiator efficiency, as in
        :func:`free_radical_network`.
    I, M : float or array-like of float
        Initiator and monomer concentrations.

    Returns
    -------
    float or ndarray

    Examples
    --------
    Doubling the initiator concentration increases the rate by only
    :math:`\sqrt2`, not 2 -- the square-root dependence:

    >>> kd, f, kp, kt, M = 1.0e-5, 0.5, 1.0e3, 1.0e7, 5.0
    >>> Rp1 = steady_state_rate_of_polymerization(kd, f, kp, kt, I=0.01, M=M)
    >>> Rp2 = steady_state_rate_of_polymerization(kd, f, kp, kt, I=0.02, M=M)
    >>> round(float(Rp2 / Rp1), 6)
    1.414214
    """
    M = np.asarray(M, dtype=np.float64)
    R_ss = steady_state_radical_concentration(kd, f, I, kt)
    result = kp * M * R_ss
    return float(result) if np.asarray(result).ndim == 0 else result


def kinetic_chain_length(kd: float, f: float, kp: float, kt: float, I, M):
    r"""Kinetic chain length :math:`\nu=R_p/R_i`: monomers consumed per radical generated.

    .. math::

        \nu = \frac{R_p}{R_i} = \frac{k_p[M]}{2\sqrt{fk_dk_t[I]}}

    directly related to the number-average degree of polymerization of
    the resulting dead polymer: :math:`\bar X_n=\nu` for termination by
    disproportionation (each radical becomes one dead chain) or
    :math:`\bar X_n=2\nu` for termination by combination (two radicals
    fuse into one dead chain) -- see Odian, *Principles of Polymerization*,
    4th ed., Ch. 3.3.2.

    Parameters
    ----------
    kd, f, kp, kt : float
        Rate constants and initiator efficiency, as in
        :func:`free_radical_network`.
    I, M : float or array-like of float
        Initiator and monomer concentrations.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> nu = kinetic_chain_length(kd=1.0e-5, f=0.5, kp=1.0e3, kt=1.0e7, I=0.01, M=5.0)
    >>> Ri = 2.0 * 0.5 * 1.0e-5 * 0.01
    >>> Rp = steady_state_rate_of_polymerization(kd=1.0e-5, f=0.5, kp=1.0e3, kt=1.0e7, I=0.01, M=5.0)
    >>> round(float(nu), 6) == round(float(Rp / Ri), 6)
    True
    """
    I = np.asarray(I, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    result = kp * M / (2.0 * np.sqrt(f * kd * kt * I))
    return float(result) if result.ndim == 0 else result
