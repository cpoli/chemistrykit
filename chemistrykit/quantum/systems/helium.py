r"""The variational helium atom: a screened effective nuclear charge for two-electron atoms and ions.

For a two-electron atom or ion of nuclear charge `Z` (H-, He, Li+, ...),
the trial wavefunction is a product of two hydrogen-like 1s orbitals
that share a *variational* effective charge :math:`\zeta`,

.. math::

    \Psi(\mathbf r_1,\mathbf r_2) = \frac{\zeta^3}{\pi a_0^3}\,e^{-\zeta(r_1+r_2)/a_0}

This gives the closed-form energy expectation value (in hartrees,
:math:`E_h=m_ee^4/(4\pi\varepsilon_0\hbar)^2`, infinite nuclear mass)

.. math::

    E(\zeta) = \left(\zeta^2 - 2Z\zeta + \tfrac{5}{8}\zeta\right)E_h

where :math:`\zeta^2` is the two electrons' kinetic energy,
:math:`-2Z\zeta` their attraction to the nucleus, and :math:`5\zeta/8` their
mutual Coulomb repulsion. Minimizing over :math:`\zeta` gives
:math:`\zeta^*=Z-5/16` and :math:`E^*=-(Z-5/16)^2E_h`: each electron
partly *screens* the nucleus from the other, by the screening constant
:math:`5/16`. For helium this is :math:`-2.8477\,E_h`, within 2% of the exact
nonrelativistic :math:`-2.9037\,E_h` (G. W. Kellner, *Z. Phys.* 44, 91
(1927); E. A. Hylleraas, *Z. Phys.* 54, 347 (1929); Levine, *Quantum
Chemistry*, 7th ed., Ch. 8.2 (variational treatment of the helium ground state)).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY

__all__ = ["HARTREE_ENERGY", "HeliumVariationalResult", "helium_like_variational_energy", "optimize_helium_like_effective_charge"]

#: The hartree, :math:`E_h=m_ee^4/((4\pi\varepsilon_0)^2\hbar^2)`, in J
#: (infinite-nuclear-mass atomic unit of energy, about 27.211 eV).
HARTREE_ENERGY = ELECTRON_MASS * ELEMENTARY_CHARGE**4 / ((4.0 * np.pi * VACUUM_PERMITTIVITY) ** 2 * HBAR**2)


def helium_like_variational_energy(zeta, Z: float = 2.0):
    r"""Energy of the screened-1s-product trial wavefunction for a two-electron atom or ion.

    :math:`E(\zeta)=(\zeta^2-2Z\zeta+5\zeta/8)\,E_h`.

    Parameters
    ----------
    zeta : float or array-like of float
        Variational effective nuclear charge (dimensionless, positive).
    Z : float, default 2.0
        True nuclear charge (2 for helium, 1 for H-, 3 for Li+, ...).

    Returns
    -------
    float or ndarray
        Energy, in J.

    Examples
    --------
    Setting :math:`\zeta=Z` recovers first-order perturbation theory,
    :math:`E=(-Z^2+5Z/8)E_h=-2.75\,E_h` for helium:

    >>> round(float(helium_like_variational_energy(2.0, Z=2.0) / HARTREE_ENERGY), 6)
    -2.75
    """
    zeta = np.asarray(zeta, dtype=np.float64)
    if np.any(zeta <= 0):
        raise ValueError("zeta must be positive")
    energy = (zeta**2 - 2.0 * Z * zeta + 0.625 * zeta) * HARTREE_ENERGY
    return float(energy) if energy.ndim == 0 else energy


@dataclass
class HeliumVariationalResult:
    """Result of :func:`optimize_helium_like_effective_charge`."""

    Z: float
    """float: True nuclear charge."""

    effective_charge: float
    r"""float: Optimal variational charge :math:`\zeta^*=Z-5/16`."""

    energy: float
    """float: Optimized variational energy, in J."""

    first_order_energy: float
    r"""float: First-order perturbation-theory energy (:math:`\zeta=Z`), in J."""

    independent_electron_energy: float
    r"""float: Energy with electron-electron repulsion ignored, :math:`-Z^2E_h`, in J."""

    @property
    def screening_constant(self) -> float:
        r"""float: :math:`Z-\zeta^*` (always :math:`5/16`)."""
        return self.Z - self.effective_charge

    @property
    def ionization_energy(self) -> float:
        r"""float: Energy to remove one electron, :math:`-Z^2E_h/2-E^*`, in J (negative means predicted unbound)."""
        return -0.5 * self.Z**2 * HARTREE_ENERGY - self.energy


def optimize_helium_like_effective_charge(Z: float = 2.0) -> HeliumVariationalResult:
    r"""Minimize :func:`helium_like_variational_energy` over the effective charge (closed form).

    :math:`dE/d\zeta=0` gives :math:`\zeta^*=Z-5/16` and
    :math:`E^*=-(Z-5/16)^2E_h`.

    Parameters
    ----------
    Z : float, default 2.0
        True nuclear charge; must exceed 5/16.

    Returns
    -------
    HeliumVariationalResult

    Examples
    --------
    Helium: :math:`\zeta^*=27/16` and :math:`E^*=-(27/16)^2E_h\approx-2.8477\,E_h`:

    >>> result = optimize_helium_like_effective_charge(2)
    >>> result.effective_charge
    1.6875
    >>> round(result.energy / HARTREE_ENERGY, 5)
    -2.84766
    """
    if Z <= 0.3125:
        raise ValueError("Z must exceed 5/16")
    zeta = float(Z) - 0.3125
    return HeliumVariationalResult(
        Z=float(Z),
        effective_charge=zeta,
        energy=-(zeta**2) * HARTREE_ENERGY,
        first_order_energy=helium_like_variational_energy(float(Z), Z=Z),
        independent_electron_energy=-(float(Z) ** 2) * HARTREE_ENERGY,
    )
