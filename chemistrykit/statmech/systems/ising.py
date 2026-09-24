r"""Exact results for the nearest-neighbor Ising model.

The Ising model places a spin :math:`s_i=\pm1` on every lattice site,
with energy :math:`E=-J\sum_{\langle ij\rangle}s_is_j-h\sum_is_i`
(coupling `J`, field `h`, both in J). It is also the lattice-gas model of
an interacting fluid or adsorbed layer (occupied/empty site
:math:`\leftrightarrow` up/down spin), which is why it sits in molecular
statistical mechanics. This module collects its three classic exact
results:

- :class:`Ising1D` -- E. Ising, "Beitrag zur Theorie des
  Ferromagnetismus," *Z. Phys.* 31, 253-258 (1925), solved here by the
  transfer matrix: no phase transition at any :math:`T>0`.
- :func:`kramers_wannier_dual_coupling` and
  :func:`ising_2d_critical_temperature` -- H. A. Kramers and G. H.
  Wannier, *Phys. Rev.* 60, 252-262 (1941): the high/low-temperature
  duality of the square lattice, which fixes :math:`T_c` exactly.
- :class:`Ising2DOnsager` -- L. Onsager, *Phys. Rev.* 65, 117-149 (1944)
  (free energy, energy, heat capacity), with C. N. Yang's spontaneous
  magnetization, *Phys. Rev.* 85, 808-816 (1952).

Formulas follow K. Huang, *Statistical Mechanics*, 2nd ed., Ch. 14-15.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipe, ellipk

from chemistrykit.constants import K_B

__all__ = ["Ising1D", "Ising2DOnsager", "kramers_wannier_dual_coupling", "ising_2d_critical_temperature"]


class Ising1D:
    r"""The one-dimensional Ising chain, solved exactly by the transfer matrix.

    With :math:`K=J/k_BT` and :math:`b=h/k_BT`, the transfer matrix
    :math:`T_{ss'}=e^{Kss'+b(s+s')/2}` has eigenvalues
    :math:`\lambda_\pm=e^K[\cosh b\pm\sqrt{\sinh^2b+e^{-4K}}]`, so a
    periodic chain of `N` spins has :math:`Z_N=\lambda_+^N+\lambda_-^N`.

    Parameters
    ----------
    coupling : float
        Exchange coupling `J`, in J (positive = ferromagnetic).
    field : float, default 0.0
        External field `h`, in J.

    Examples
    --------
    In zero field the magnetization vanishes at every temperature -- the
    one-dimensional chain never orders:

    >>> chain = Ising1D(coupling=1.380649e-23 * 100.0)
    >>> float(chain.magnetization(50.0))
    0.0

    but the correlation length grows as :math:`\xi\approx e^{2K}/2` on cooling:

    >>> round(float(chain.correlation_length(50.0)), 2)
    27.3
    """

    def __init__(self, coupling: float, field: float = 0.0):
        self.coupling = float(coupling)
        self.field = float(field)

    def _K_b(self, T):
        T = np.asarray(T, dtype=float)
        if np.any(T <= 0):
            raise ValueError("T must be positive")
        return self.coupling / (K_B * T), self.field / (K_B * T)

    def transfer_matrix_eigenvalues(self, T):
        r"""Return :math:`(\lambda_+,\lambda_-)` of the 2x2 transfer matrix at temperature `T`.

        Parameters
        ----------
        T : float
            Absolute temperature, in K.

        Returns
        -------
        tuple of float
        """
        K, b = self._K_b(T)
        root = np.sqrt(np.sinh(b) ** 2 + np.exp(-4.0 * K))
        return np.exp(K) * (np.cosh(b) + root), np.exp(K) * (np.cosh(b) - root)

    def partition_function(self, T, N: int):
        r"""Return :math:`Z_N=\lambda_+^N+\lambda_-^N` for a periodic ring of `N` spins.

        Parameters
        ----------
        T : float
            Absolute temperature, in K.
        N : int
            Number of spins.

        Returns
        -------
        float
        """
        lam_plus, lam_minus = self.transfer_matrix_eigenvalues(T)
        return lam_plus**N + lam_minus**N

    def free_energy_per_spin(self, T):
        r"""Return the thermodynamic-limit free energy per spin, :math:`f=-k_BT\ln\lambda_+`, in J."""
        lam_plus, _ = self.transfer_matrix_eigenvalues(T)
        return -K_B * np.asarray(T, dtype=float) * np.log(lam_plus)

    def magnetization(self, T):
        r"""Return the magnetization per spin, :math:`m=\sinh b/\sqrt{\sinh^2b+e^{-4K}}`."""
        K, b = self._K_b(T)
        return np.sinh(b) / np.sqrt(np.sinh(b) ** 2 + np.exp(-4.0 * K))

    def correlation_length(self, T):
        r"""Return the zero-field correlation length in lattice spacings, :math:`\xi=-1/\ln\tanh K`.

        Uses the zero-field ratio :math:`\lambda_-/\lambda_+=\tanh K`
        (the stored `field` is ignored).
        """
        K, _ = self._K_b(T)
        return -1.0 / np.log(np.tanh(K))

    def heat_capacity_per_spin(self, T):
        r"""Return the zero-field heat capacity per spin, :math:`c=k_BK^2\operatorname{sech}^2K`, in J/K."""
        K, _ = self._K_b(T)
        return K_B * K**2 / np.cosh(K) ** 2


def kramers_wannier_dual_coupling(K):
    r"""Return the Kramers-Wannier dual coupling :math:`K^*=-\tfrac12\ln\tanh K`.

    The square-lattice Ising partition function at reduced coupling
    :math:`K=J/k_BT` equals (up to an analytic prefactor) the one at
    :math:`K^*`, equivalently :math:`\sinh 2K\,\sinh 2K^*=1`: high
    temperature maps onto low temperature.

    Parameters
    ----------
    K : float or array_like
        Reduced coupling :math:`J/k_BT` (positive).

    Returns
    -------
    float or numpy.ndarray

    Examples
    --------
    The map is an involution, and its fixed point is the critical coupling
    :math:`K_c=\tfrac12\ln(1+\sqrt2)`:

    >>> round(float(kramers_wannier_dual_coupling(kramers_wannier_dual_coupling(0.3))), 12)
    0.3
    >>> Kc = 0.5 * np.log(1.0 + np.sqrt(2.0))
    >>> round(float(kramers_wannier_dual_coupling(Kc) - Kc), 12)
    0.0
    """
    return -0.5 * np.log(np.tanh(np.asarray(K, dtype=float)))


def ising_2d_critical_temperature(coupling: float) -> float:
    r"""Return the exact square-lattice critical temperature :math:`T_c=2J/[k_B\ln(1+\sqrt2)]`.

    Parameters
    ----------
    coupling : float
        Exchange coupling `J`, in J.

    Returns
    -------
    float
        Critical temperature, in K.

    Examples
    --------
    In units of :math:`J/k_B`, :math:`T_c\approx2.269`:

    >>> round(ising_2d_critical_temperature(1.380649e-23), 4)
    2.2692
    """
    return float(2.0 * coupling / (K_B * np.log(1.0 + np.sqrt(2.0))))


class Ising2DOnsager:
    r"""Onsager's exact solution of the zero-field square-lattice Ising model.

    With :math:`K=J/k_BT` and :math:`\kappa=2\sinh2K/\cosh^22K`:

    .. math::

        -\beta f = \ln(2\cosh2K)+\frac1\pi\int_0^{\pi/2}
        \ln\tfrac12\left(1+\sqrt{1-\kappa^2\sin^2\phi}\right)d\phi

    and, for :math:`T<T_c`, Yang's spontaneous magnetization
    :math:`m=[1-\sinh^{-4}2K]^{1/8}`.

    Parameters
    ----------
    coupling : float
        Exchange coupling `J`, in J (positive).

    Examples
    --------
    >>> model = Ising2DOnsager(coupling=1.380649e-23)  # J/k_B = 1 K
    >>> round(model.critical_temperature, 4)
    2.2692
    >>> round(float(model.spontaneous_magnetization(1.0)), 4)
    0.9993
    >>> float(model.spontaneous_magnetization(3.0))
    0.0
    >>> round(float(model.internal_energy_per_spin(model.critical_temperature) / model.coupling), 6)
    -1.414214
    """

    def __init__(self, coupling: float):
        if coupling <= 0:
            raise ValueError("coupling must be positive")
        self.coupling = float(coupling)

    @property
    def critical_temperature(self) -> float:
        """The exact critical temperature, in K. See :func:`ising_2d_critical_temperature`."""
        return ising_2d_critical_temperature(self.coupling)

    def _K(self, T):
        T = np.asarray(T, dtype=float)
        if np.any(T <= 0):
            raise ValueError("T must be positive")
        return self.coupling / (K_B * T)

    def free_energy_per_spin(self, T):
        """Return the free energy per spin, in J (Onsager's single-integral form)."""

        def one(t):
            K = self.coupling / (K_B * t)
            kappa = 2.0 * np.sinh(2.0 * K) / np.cosh(2.0 * K) ** 2
            integral, _ = quad(lambda phi: np.log(0.5 * (1.0 + np.sqrt(max(0.0, 1.0 - kappa**2 * np.sin(phi) ** 2)))), 0.0, np.pi / 2, limit=200)
            return -K_B * t * (np.log(2.0 * np.cosh(2.0 * K)) + integral / np.pi)

        self._K(T)
        values = np.vectorize(one)(np.asarray(T, dtype=float))
        return float(values) if values.ndim == 0 else values

    def internal_energy_per_spin(self, T):
        r"""Return :math:`u=-J\coth2K\,[1+\tfrac2\pi(2\tanh^22K-1)K_1(\kappa)]`, in J.

        :math:`K_1` is the complete elliptic integral of the first kind; its
        logarithmic divergence at :math:`T_c` is multiplied by zero, so `u`
        is continuous there (:math:`u(T_c)=-\sqrt2J`).
        """
        K = self._K(T)
        kappa = 2.0 * np.sinh(2.0 * K) / np.cosh(2.0 * K) ** 2
        kappa_p = 2.0 * np.tanh(2.0 * K) ** 2 - 1.0
        with np.errstate(invalid="ignore", divide="ignore"):
            term = kappa_p * ellipk(np.minimum(kappa**2, 1.0))
        term = np.where(np.isfinite(term), term, 0.0)
        return -self.coupling / np.tanh(2.0 * K) * (1.0 + 2.0 / np.pi * term)

    def heat_capacity_per_spin(self, T):
        r"""Return the heat capacity per spin, in J/K (diverges logarithmically at :math:`T_c`).

        .. math::

            \frac{c}{k_B}=\frac2\pi(K\coth2K)^2\Big\{2K_1(\kappa)-2E_1(\kappa)
            -(1-\kappa')\big[\tfrac\pi2+\kappa'K_1(\kappa)\big]\Big\}

        with :math:`\kappa'=2\tanh^22K-1` and :math:`E_1` the complete
        elliptic integral of the second kind.
        """
        K = self._K(T)
        kappa = 2.0 * np.sinh(2.0 * K) / np.cosh(2.0 * K) ** 2
        kappa_p = 2.0 * np.tanh(2.0 * K) ** 2 - 1.0
        m = np.minimum(kappa**2, 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            K1, E1 = ellipk(m), ellipe(m)
            bracket = 2.0 * K1 - 2.0 * E1 - (1.0 - kappa_p) * (np.pi / 2.0 + kappa_p * K1)
        bracket = np.where(m >= 1.0, np.inf, bracket)
        return K_B * 2.0 / np.pi * (K / np.tanh(2.0 * K)) ** 2 * bracket

    def spontaneous_magnetization(self, T):
        r"""Return Yang's spontaneous magnetization per spin, :math:`[1-\sinh^{-4}2K]^{1/8}` below :math:`T_c`, else 0."""
        K = self._K(T)
        with np.errstate(divide="ignore", invalid="ignore"):
            inner = 1.0 - np.sinh(2.0 * K) ** -4.0
        return np.where(inner > 0.0, np.abs(inner) ** 0.125, 0.0)
