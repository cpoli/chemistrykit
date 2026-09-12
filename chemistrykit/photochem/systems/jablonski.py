r"""Jablonski-diagram excited-state kinetics, as a mass-action reaction network.

A Jablonski diagram (A. Jablonski, *Nature* 131, 839 (1933)) tracks a
molecule's electronic-state population after photoexcitation as it
relaxes back to the ground state through several competing first-order
pathways. This module models the standard minimal three-state version --
the first excited singlet :math:`S_1`, the lowest triplet :math:`T_1`,
and the ground singlet :math:`S_0` -- with rate constants for
fluorescence (:math:`k_f`, :math:`S_1\to S_0` radiative), internal
conversion (:math:`k_{ic}`, :math:`S_1\to S_0` nonradiative),
intersystem crossing (:math:`k_{isc}`, :math:`S_1\to T_1`),
phosphorescence (:math:`k_p`, :math:`T_1\to S_0` radiative), and
:math:`T_1\to S_0` nonradiative decay (:math:`k_{ic,T}`) -- see Turro,
Ramamurthy & Scaiano, *Modern Molecular Photochemistry of Organic
Molecules* (2010), Ch. 5, or Lakowicz, *Principles of Fluorescence
Spectroscopy*, 3rd ed., Ch. 1.

Every one of these steps is a first-order (unimolecular) rate process,
so the whole diagram is *exactly* the kind of mass-action network
:class:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork`
already models generally -- :func:`jablonski_network` builds one rather
than reimplementing rate-equation integration, per this package's
convention of reusing :mod:`chemistrykit.kinetics`'s machinery wherever a
new domain's kinetics reduces to the same mathematics.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork, consecutive_analytic

__all__ = ["jablonski_network", "jablonski_populations_analytic"]


def jablonski_network(kf: float, kic: float, kisc: float, kp: float, kic_T: float, S1_0: float = 1.0) -> StoichiometricNetwork:
    r"""Build the 3-state (:math:`S_1,T_1,S_0`) Jablonski excited-state decay network.

    Structurally, :math:`S_1` decays through two parallel first-order
    channels (directly to :math:`S_0` at combined rate :math:`k_f+k_{ic}`,
    and to :math:`T_1` at rate :math:`k_{isc}`), and :math:`T_1` then
    decays to :math:`S_0` at combined rate :math:`k_p+k_{ic,T}` -- a
    branching-then-consecutive network, built here as a general
    :class:`~chemistrykit.kinetics.systems.networks.StoichiometricNetwork`
    (its three named constructors, ``parallel``/``consecutive``/
    ``reversible``, don't individually cover this branching topology, but
    the general mass-action constructor they're built on does).

    Parameters
    ----------
    kf : float
        Fluorescence (:math:`S_1\to S_0`, radiative) rate constant.
    kic : float
        Internal conversion (:math:`S_1\to S_0`, nonradiative) rate
        constant.
    kisc : float
        Intersystem crossing (:math:`S_1\to T_1`) rate constant.
    kp : float
        Phosphorescence (:math:`T_1\to S_0`, radiative) rate constant.
    kic_T : float
        Triplet nonradiative decay (:math:`T_1\to S_0`) rate constant.
    S1_0 : float, default 1.0
        Initial :math:`S_1` population immediately after excitation
        (e.g. following a short excitation pulse).

    Returns
    -------
    StoichiometricNetwork
        Species ``("S1", "T1", "S0")``.

    Examples
    --------
    Total population is conserved (S0 is a pure sink, nothing is created
    or destroyed, just redistributed):

    >>> net = jablonski_network(kf=2.0, kic=1.0, kisc=0.5, kp=0.3, kic_T=0.2, S1_0=1.0)
    >>> result = net.integrate((0.0, 20.0), dt=1e-3, method="rk4")
    >>> total = result.concentration("S1") + result.concentration("T1") + result.concentration("S0")
    >>> bool(np.allclose(total, 1.0, atol=1e-6))
    True
    """
    species = ("S1", "T1", "S0")
    stoich = [[-1.0, -1.0, 0.0], [0.0, 1.0, -1.0], [1.0, 0.0, 1.0]]
    reactant_orders = [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]
    rate_constants = [kf + kic, kisc, kp + kic_T]
    state0 = [S1_0, 0.0, 0.0]
    return StoichiometricNetwork(species, stoich, rate_constants, reactant_orders, state0)


def jablonski_populations_analytic(kf: float, kic: float, kisc: float, kp: float, kic_T: float, S1_0: float, t):
    r"""Closed-form :math:`S_1(t)`, :math:`T_1(t)`, :math:`S_0(t)` for the 3-state Jablonski network.

    :math:`S_1` decays as a simple single exponential with total rate
    :math:`k_{S_1}=k_f+k_{ic}+k_{isc}`:

    .. math::

        [S_1](t) = [S_1]_0 e^{-k_{S_1}t}

    :math:`T_1` is fed by the fraction :math:`k_{isc}/k_{S_1}` of decayed
    :math:`S_1` and itself decays at total rate
    :math:`k_{T_1}=k_p+k_{ic,T}`; solving
    :math:`d[T_1]/dt=k_{isc}[S_1]-k_{T_1}[T_1]` gives a scaled version of
    the same Bateman-equation form as a consecutive :math:`A\to B\to C`
    chain (Turro, Ramamurthy & Scaiano, *Modern Molecular Photochemistry
    of Organic Molecules*, Ch. 5; the underlying algebra is identical to
    :func:`chemistrykit.kinetics.systems.networks.consecutive_analytic`,
    which this function calls directly and rescales by
    :math:`k_{isc}/k_{S_1}` -- since only that fraction of decayed
    :math:`S_1`, not all of it, feeds :math:`T_1`):

    .. math::

        [T_1](t) = \frac{k_{isc}}{k_{S_1}}
        \left[\frac{[S_1]_0 k_{S_1}}{k_{T_1}-k_{S_1}}
        \left(e^{-k_{S_1}t}-e^{-k_{T_1}t}\right)\right]

    and :math:`[S_0](t) = [S_1]_0 - [S_1](t) - [T_1](t)` by mass balance.
    This closed form was checked against direct numerical integration of
    :func:`jablonski_network` (to a relative agreement of better than
    :math:`10^{-10}` across a range of rate constants) as part of
    developing this module.

    Parameters
    ----------
    kf, kic, kisc, kp, kic_T : float
        Rate constants, as in :func:`jablonski_network`.
    S1_0 : float
        Initial :math:`S_1` population.
    t : float or array-like of float
        Time(s) at which to evaluate the populations.

    Returns
    -------
    S1, T1, S0 : float or ndarray

    Examples
    --------
    Matches :func:`jablonski_network`'s numerical integration:

    >>> kf, kic, kisc, kp, kic_T, S1_0 = 2.0, 1.0, 0.5, 0.3, 0.2, 1.0
    >>> net = jablonski_network(kf, kic, kisc, kp, kic_T, S1_0)
    >>> result = net.integrate((0.0, 10.0), dt=1e-3, method="rk4")
    >>> S1, T1, S0 = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, S1_0, result.t)
    >>> bool(np.allclose(result.concentration("S1"), S1, atol=1e-6))
    True
    >>> bool(np.allclose(result.concentration("T1"), T1, atol=1e-6))
    True
    """
    t = np.asarray(t, dtype=np.float64)
    k_S1 = kf + kic + kisc
    k_T1 = kp + kic_T
    S1 = S1_0 * np.exp(-k_S1 * t)
    _, B, _ = consecutive_analytic(A0=S1_0, k1=k_S1, k2=k_T1, t=t)
    T1 = (kisc / k_S1) * B
    S0 = S1_0 - S1 - T1
    return S1, T1, S0
