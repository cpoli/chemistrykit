r"""The photochemical H2 + Cl2 chain reaction, and its super-unity quantum yield.

Bodenstein's 1913 measurement that each photon absorbed by a
hydrogen/chlorine mixture forms of order :math:`10^5`-:math:`10^6` HCl
molecules, and Nernst's 1918 atom-chain explanation of it, are modeled
here with the minimal Nernst mechanism (Laidler, *Chemical Kinetics*,
3rd ed. (1987), Ch. 8; Atkins & de Paula, *Physical Chemistry*, Ch. on
chain reactions):

.. math::

    \mathrm{Cl_2} + h\nu &\xrightarrow{j} 2\,\mathrm{Cl}
    \qquad \text{(initiation, } I_{abs}=j[\mathrm{Cl_2}]\text{)}\\
    \mathrm{Cl} + \mathrm{H_2} &\xrightarrow{k_2} \mathrm{HCl} + \mathrm{H}
    \qquad \text{(propagation)}\\
    \mathrm{H} + \mathrm{Cl_2} &\xrightarrow{k_3} \mathrm{HCl} + \mathrm{Cl}
    \qquad \text{(propagation)}\\
    \mathrm{Cl} + \mathrm{Cl} &\xrightarrow{k_t} \mathrm{Cl_2}
    \qquad \text{(termination)}

Applying the steady-state approximation to both atoms gives
:math:`[\mathrm{Cl}]_{ss}=\sqrt{I_{abs}/k_t}` and a rate of HCl formation
:math:`2k_2[\mathrm{Cl}][\mathrm{H_2}]`, so the quantum yield (HCl
molecules per absorbed photon) is

.. math::

    \Phi_{HCl} = \frac{2k_2[\mathrm{H_2}]}{\sqrt{k_t I_{abs}}},

which is unbounded above: long chains (fast propagation, slow
termination, weak light) give :math:`\Phi\gg1`, the historically real
exception to the Stark-Einstein bound :math:`\Phi\le1`.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork

__all__ = ["hydrogen_chlorine_chain_network", "chain_quantum_yield"]


def hydrogen_chlorine_chain_network(j: float, k2: float, k3: float, kt: float, H2_0: float = 1.0, Cl2_0: float = 1.0) -> StoichiometricNetwork:
    r"""Build Nernst's photochemical H2/Cl2 atom-chain mechanism as a mass-action network.

    Parameters
    ----------
    j : float
        Photolysis rate constant of :math:`\mathrm{Cl_2}` (so the
        absorbed photon rate per volume is :math:`I_{abs}=j[\mathrm{Cl_2}]`,
        assuming each absorbed photon dissociates one :math:`\mathrm{Cl_2}`).
    k2 : float
        :math:`\mathrm{Cl}+\mathrm{H_2}\to\mathrm{HCl}+\mathrm{H}` rate constant.
    k3 : float
        :math:`\mathrm{H}+\mathrm{Cl_2}\to\mathrm{HCl}+\mathrm{Cl}` rate constant.
    kt : float
        :math:`\mathrm{Cl}+\mathrm{Cl}\to\mathrm{Cl_2}` termination rate
        constant (mass-action rate :math:`k_t[\mathrm{Cl}]^2`).
    H2_0, Cl2_0 : float, default 1.0
        Initial concentrations.

    Returns
    -------
    StoichiometricNetwork
        Species ``("H2", "Cl2", "Cl", "H", "HCl")``.

    Examples
    --------
    Hydrogen atoms are conserved between H2, H and HCl:

    >>> net = hydrogen_chlorine_chain_network(j=1e-4, k2=10.0, k3=100.0, kt=100.0)
    >>> res = net.integrate((0.0, 5.0), dt=1e-3, method="rk4")
    >>> H_total = 2 * res.concentration("H2") + res.concentration("H") + res.concentration("HCl")
    >>> bool(np.allclose(H_total, 2.0))
    True
    """
    species = ("H2", "Cl2", "Cl", "H", "HCl")
    # columns: initiation, Cl + H2, H + Cl2, Cl + Cl
    stoich = [
        [0.0, -1.0, 0.0, 0.0],
        [-1.0, 0.0, -1.0, 1.0],
        [2.0, -1.0, 1.0, -2.0],
        [0.0, 1.0, -1.0, 0.0],
        [0.0, 1.0, 1.0, 0.0],
    ]
    orders = [
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 2.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    state0 = [H2_0, Cl2_0, 0.0, 0.0, 0.0]
    return StoichiometricNetwork(species, stoich, [j, k2, k3, kt], orders, state0)


def chain_quantum_yield(k2: float, kt: float, H2, I_abs):
    r"""Steady-state quantum yield of the H2/Cl2 chain: :math:`\Phi = 2k_2[\mathrm{H_2}]/\sqrt{k_t I_{abs}}`.

    Parameters
    ----------
    k2 : float
        :math:`\mathrm{Cl}+\mathrm{H_2}` propagation rate constant.
    kt : float
        :math:`\mathrm{Cl}+\mathrm{Cl}` termination rate constant.
    H2 : float or array-like of float
        Hydrogen concentration.
    I_abs : float or array-like of float
        Absorbed photon rate per unit volume (same concentration/time
        units as the rate constants).

    Returns
    -------
    float or ndarray
        HCl molecules formed per photon absorbed -- not bounded by 1.

    Examples
    --------
    Halving the light intensity by a factor of 4 doubles the quantum
    yield (the chain lengthens as :math:`I_{abs}^{-1/2}`):

    >>> phi1 = chain_quantum_yield(k2=10.0, kt=100.0, H2=1.0, I_abs=1e-4)
    >>> phi2 = chain_quantum_yield(k2=10.0, kt=100.0, H2=1.0, I_abs=0.25e-4)
    >>> round(phi1, 6), round(phi2 / phi1, 6)
    (200.0, 2.0)
    """
    result = 2.0 * k2 * np.asarray(H2, dtype=np.float64) / np.sqrt(kt * np.asarray(I_abs, dtype=np.float64))
    return float(result) if result.ndim == 0 else result
