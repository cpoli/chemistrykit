r"""Hydrogen-like (one-electron) atoms: exact energies, radial wavefunctions, and orbital shapes.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 8.1-8.2, or
Levine, *Quantum Chemistry*, 7th ed., Ch. 6.6-6.7, for the exact analytic
solution of the Schrodinger equation for one electron in the Coulomb field
of a nucleus of charge `Z`. The angular part (spherical harmonics) and the
radial part's associated Laguerre polynomials are taken from
:mod:`scipy.special` rather than re-derived, per the recurring pattern
elsewhere in this subpackage of reusing well-tested special-function
implementations for the purely mathematical piece of an exact solution.
"""

from __future__ import annotations

import numpy as np
import scipy.special as sp
from scipy.integrate import quad

from chemistrykit.constants import ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY
from chemistrykit.quantum.core.base_system import QuantumSystem

__all__ = ["HydrogenLikeAtom"]


def _spherical_harmonic(l: int, m: int, theta, phi):
    """Return :math:`Y_l^m(\\theta,\\phi)` (`theta` polar, `phi` azimuthal), across scipy versions.

    scipy 1.15 renamed/replaced ``scipy.special.sph_harm`` (removed in
    1.18) with ``sph_harm_y``, which also swapped nothing about the
    theta/phi *meaning* but reordered the call signature and argument
    names; this shim keeps :class:`HydrogenLikeAtom` working across the
    ``scipy>=1.10`` range this package supports.
    """
    if hasattr(sp, "sph_harm_y"):
        return sp.sph_harm_y(l, m, theta, phi)
    return sp.sph_harm(m, l, phi, theta)  # pragma: no cover - only exercised on scipy < 1.15


class HydrogenLikeAtom(QuantumSystem):
    r"""A one-electron atom/ion of nuclear charge `Z` (H, He+, Li2+, ...).

    .. math::

        E_n = -\frac{Z^2\mu e^4}{2(4\pi\varepsilon_0)^2\hbar^2n^2}, \qquad n=1,2,3,\dots

    (Atkins & de Paula, *Physical Chemistry*, 11th ed., eq. 8.2.) The
    radial wavefunctions are

    .. math::

        R_{n,l}(r) = \sqrt{\left(\frac{2Z}{na_0}\right)^3\frac{(n-l-1)!}{2n\,(n+l)!}}
                     e^{-\rho/2}\rho^lL_{n-l-1}^{2l+1}(\rho), \qquad \rho=\frac{2Zr}{na_0}

    with :math:`a_0=4\pi\varepsilon_0\hbar^2/(\mu e^2)` the (reduced-mass)
    Bohr radius and :math:`L_{n-l-1}^{2l+1}` an associated Laguerre
    polynomial (Atkins & de Paula, *Physical Chemistry*, 11th ed., Table
    8.1; Levine, *Quantum Chemistry*, 7th ed., eq. 6.94), and the full
    wavefunction is :math:`\psi_{n,l,m}(r,\theta,\phi)=R_{n,l}(r)Y_l^m(\theta,\phi)`.

    By default `reduced_mass` is the bare electron mass -- the
    infinite-nuclear-mass approximation, giving the textbook
    :math:`-13.6\,\text{eV}` hydrogen ground state to 3 significant
    figures. Passing the true two-body reduced mass
    :math:`\mu=m_em_{nuc}/(m_e+m_{nuc})` instead removes this
    approximation and reproduces hydrogen's ground-state energy to
    5-6 figures (:math:`-13.598\,\text{eV}`) -- the small residual
    difference is exactly the nuclear recoil the infinite-mass
    approximation neglects.

    Parameters
    ----------
    Z : int, default 1
        Nuclear charge (in units of the elementary charge).
    reduced_mass : float, default :data:`chemistrykit.constants.ELECTRON_MASS`
        Reduced mass of the electron-nucleus system, in kg.

    Examples
    --------
    Hydrogen's ground-state ionization energy is the textbook 13.6 eV:

    >>> h_atom = HydrogenLikeAtom(Z=1)
    >>> round(float(-h_atom.energy(1) / 1.602176634e-19), 1)
    13.6

    He+ (Z=2) is bound 4x more tightly than H at the same `n` (:math:`E\propto Z^2`):

    >>> he_plus = HydrogenLikeAtom(Z=2)
    >>> round(float(he_plus.energy(1) / h_atom.energy(1)), 6)
    4.0
    """

    def __init__(self, Z: int = 1, reduced_mass: float = ELECTRON_MASS):
        if Z < 1:
            raise ValueError("Z must be >= 1")
        if reduced_mass <= 0:
            raise ValueError("reduced_mass must be positive")
        self.Z = int(Z)
        self.reduced_mass = float(reduced_mass)

    @property
    def bohr_radius(self) -> float:
        r"""float: The (reduced-mass) Bohr radius :math:`a_0=4\pi\varepsilon_0\hbar^2/(\mu e^2)`, in m."""
        return 4.0 * np.pi * VACUUM_PERMITTIVITY * HBAR**2 / (self.reduced_mass * ELEMENTARY_CHARGE**2)

    def energy(self, n):
        r"""Return :math:`E_n=-Z^2\mu e^4/(2(4\pi\varepsilon_0)^2\hbar^2n^2)`.

        Parameters
        ----------
        n : int or array-like of int
            Principal quantum number(s), each >= 1.

        Returns
        -------
        float or ndarray
            Energy, in J (negative -- bound).
        """
        n = np.asarray(n, dtype=np.float64)
        if np.any(n < 1):
            raise ValueError("n must be >= 1")
        numerator = self.Z**2 * self.reduced_mass * ELEMENTARY_CHARGE**4
        denominator = 2.0 * (4.0 * np.pi * VACUUM_PERMITTIVITY) ** 2 * HBAR**2 * n**2
        return -numerator / denominator

    def radial_wavefunction(self, r, n: int, l: int):
        r"""Return the radial wavefunction :math:`R_{n,l}(r)`.

        Parameters
        ----------
        r : float or array-like of float
            Radial distance(s) from the nucleus, in m.
        n : int
            Principal quantum number, >= 1.
        l : int
            Orbital angular momentum quantum number, ``0 <= l <= n - 1``.

        Returns
        -------
        float or ndarray
            Amplitude, in m^-3/2.

        Examples
        --------
        The 1s radial wavefunction is a simple decaying exponential,
        :math:`R_{1,0}(r)=2a_0^{-3/2}e^{-r/a_0}` (Atkins & de Paula,
        *Physical Chemistry*, 11th ed., Table 8.1):

        >>> h_atom = HydrogenLikeAtom(Z=1)
        >>> a0 = h_atom.bohr_radius
        >>> R10 = h_atom.radial_wavefunction(a0, n=1, l=0)
        >>> expected = 2.0 * a0 ** -1.5 * np.exp(-1.0)
        >>> round(float(R10 / expected), 6)
        1.0
        """
        if not (0 <= l <= n - 1):
            raise ValueError("l must satisfy 0 <= l <= n - 1")
        r = np.asarray(r, dtype=np.float64)
        a0 = self.bohr_radius
        rho = 2.0 * self.Z * r / (n * a0)
        norm = np.sqrt((2.0 * self.Z / (n * a0)) ** 3 * sp.factorial(n - l - 1) / (2.0 * n * sp.factorial(n + l)))
        laguerre = sp.genlaguerre(n - l - 1, 2 * l + 1)
        return norm * np.exp(-rho / 2.0) * rho**l * laguerre(rho)

    def radial_distribution_function(self, r, n: int, l: int):
        r"""Return the radial distribution function :math:`P(r)=r^2R_{n,l}(r)^2`.

        The probability of finding the electron in a thin spherical shell
        between `r` and `r+dr`, obtained by integrating :math:`|\psi|^2`
        over all angles (Atkins & de Paula, *Physical Chemistry*, 11th
        ed., eq. 8.4).

        Parameters
        ----------
        r : float or array-like of float
        n : int
        l : int

        Returns
        -------
        float or ndarray
            Probability density, in m^-1.
        """
        r = np.asarray(r, dtype=np.float64)
        return r**2 * self.radial_wavefunction(r, n, l) ** 2

    def angular_wavefunction(self, theta, phi, l: int, m: int):
        r"""Return the (generally complex) spherical harmonic :math:`Y_l^m(\theta,\phi)`.

        Parameters
        ----------
        theta : float or array-like of float
            Polar angle, in radians, in ``[0, pi]``.
        phi : float or array-like of float
            Azimuthal angle, in radians, in ``[0, 2*pi)``.
        l : int
            Orbital angular momentum quantum number, >= 0.
        m : int
            Magnetic quantum number, ``-l <= m <= l``.

        Returns
        -------
        complex or ndarray of complex
        """
        if not (-l <= m <= l):
            raise ValueError("m must satisfy -l <= m <= l")
        return _spherical_harmonic(l, m, theta, phi)

    def wavefunction(self, r, theta, phi, n: int, l: int, m: int):
        r"""Return the full (generally complex) wavefunction :math:`\psi_{n,l,m}(r,\theta,\phi)=R_{n,l}(r)Y_l^m(\theta,\phi)`.

        Parameters
        ----------
        r : float or array-like of float
            Radial distance, in m.
        theta : float or array-like of float
            Polar angle, in radians.
        phi : float or array-like of float
            Azimuthal angle, in radians.
        n, l, m : int
            Quantum numbers, ``0 <= l <= n - 1``, ``-l <= m <= l``.

        Returns
        -------
        complex or ndarray of complex
            Amplitude, in m^-3/2.
        """
        return self.radial_wavefunction(r, n, l) * self.angular_wavefunction(theta, phi, l, m)

    def check_radial_normalization(self, n: int, l: int, r_max_bohr_radii: float = 60.0) -> float:
        r"""Numerically verify :math:`\int_0^\infty R_{n,l}(r)^2r^2\,dr=1`.

        A direct, formula-agnostic sanity check of :meth:`radial_wavefunction`
        (independent of any particular normalization-constant convention):
        integrates the radial distribution function
        (:meth:`radial_distribution_function`) out to `r_max_bohr_radii`
        Bohr radii, where the (exponentially decaying) integrand is
        already negligible.

        Parameters
        ----------
        n, l : int
        r_max_bohr_radii : float, default 60.0
            Upper integration limit, in units of :attr:`bohr_radius`.

        Returns
        -------
        float
            Should equal 1.0 to within numerical-quadrature error.

        Examples
        --------
        >>> h_atom = HydrogenLikeAtom(Z=1)
        >>> round(h_atom.check_radial_normalization(n=1, l=0), 6)
        1.0
        >>> round(h_atom.check_radial_normalization(n=2, l=1), 6)
        1.0
        """
        a0 = self.bohr_radius
        integral, _ = quad(lambda r: self.radial_distribution_function(r, n, l), 0.0, r_max_bohr_radii * a0, limit=200)
        return float(integral)
