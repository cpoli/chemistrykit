r"""Photoinduced electron-transfer quenching: the Rehm-Weller equation.

D. Rehm and A. Weller, *Isr. J. Chem.* 8, 259 (1970), showed that the
rate constant for fluorescence quenching by electron transfer in
acetonitrile follows from the free energy of the electron-transfer step,

.. math::

    \Delta G_{ET} = E_{ox}(D) - E_{red}(A) - E_{00} + w,

(potentials in V, energies in eV, :math:`E_{00}` the excited-state
energy, :math:`w` the Coulombic work term), through the empirical
relation

.. math::

    k_q = \frac{k_d}{1 + \frac{k_{-d}}{k_{e}^0}\left[
    \exp(\Delta G^\ddagger/RT) + \exp(\Delta G_{ET}/RT)\right]},\qquad
    \Delta G^\ddagger = \frac{\Delta G_{ET}}{2} +
    \left[\left(\frac{\Delta G_{ET}}{2}\right)^2 +
    \left(\Delta G^\ddagger(0)\right)^2\right]^{1/2},

with :math:`k_d=2.0\times10^{10}` M\ :sup:`-1` s\ :sup:`-1`,
:math:`k_{-d}/k_e^0=0.25` and :math:`\Delta G^\ddagger(0)=2.4` kcal/mol
(0.104 eV) in their original fit (see also Turro, Ramamurthy & Scaiano,
*Modern Molecular Photochemistry of Organic Molecules*, Ch. 7). Unlike
Marcus theory, the Rehm-Weller curve plateaus at the diffusion limit for
strongly exergonic transfer instead of turning over into an inverted
region.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, K_B

__all__ = ["rehm_weller_free_energy", "rehm_weller_quenching_rate"]


def rehm_weller_free_energy(E_ox_donor: float, E_red_acceptor: float, E00: float, work_term: float = 0.0) -> float:
    r"""Free energy of photoinduced electron transfer, :math:`\Delta G_{ET}=E_{ox}(D)-E_{red}(A)-E_{00}+w` (eV).

    Parameters
    ----------
    E_ox_donor : float
        Oxidation potential of the donor, in V.
    E_red_acceptor : float
        Reduction potential of the acceptor, in V (same reference).
    E00 : float
        Excitation energy of the excited partner, in eV.
    work_term : float, default 0.0
        Coulombic ion-pair work term, in eV (small in polar solvents).

    Returns
    -------
    float
        :math:`\Delta G_{ET}`, in eV (negative means exergonic).

    Examples
    --------
    >>> round(rehm_weller_free_energy(1.2, -1.9, 3.3), 6)
    -0.2
    """
    return E_ox_donor - E_red_acceptor - E00 + work_term


def rehm_weller_quenching_rate(delta_G, k_diff: float = 2.0e10, ratio: float = 0.25, dG0_dagger: float = 0.104, T: float = 298.15):
    r"""Rehm-Weller electron-transfer quenching rate constant :math:`k_q(\Delta G_{ET})`.

    Parameters
    ----------
    delta_G : float or array-like of float
        Electron-transfer free energy, in eV.
    k_diff : float, default 2.0e10
        Diffusion rate constant :math:`k_d`, in 1/(M*s).
    ratio : float, default 0.25
        :math:`k_{-d}/k_e^0`.
    dG0_dagger : float, default 0.104
        Intrinsic barrier :math:`\Delta G^\ddagger(0)`, in eV.
    T : float, default 298.15
        Temperature, in K.

    Returns
    -------
    float or ndarray
        :math:`k_q`, in 1/(M*s).

    Examples
    --------
    Strongly exergonic transfer approaches the plateau
    :math:`k_d/(1+0.25)=1.6\times10^{10}`, while thermoneutral transfer
    is more than ten times slower:

    >>> round(rehm_weller_quenching_rate(-1.5) / 1e10, 2)
    1.5
    >>> round(rehm_weller_quenching_rate(-50.0) / 1e10, 2)
    1.6
    >>> round(rehm_weller_quenching_rate(0.0) / 1e10, 3)
    0.128
    """
    dG = np.asarray(delta_G, dtype=np.float64)
    kT_eV = K_B * T / ELEMENTARY_CHARGE
    dG_dagger = dG / 2.0 + np.sqrt((dG / 2.0) ** 2 + dG0_dagger**2)
    result = k_diff / (1.0 + ratio * (np.exp(dG_dagger / kT_eV) + np.exp(dG / kT_eV)))
    return float(result) if result.ndim == 0 else result
