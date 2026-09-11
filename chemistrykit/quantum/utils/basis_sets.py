r"""Minimal s-type Gaussian-primitive integrals, for the toy Hartree-Fock-style LCAO solver.

Only what :class:`chemistrykit.quantum.systems.hartree_fock.H2PlusVariational`
needs: the overlap, kinetic-energy, and (one-electron) nuclear-attraction
integrals between two normalized, spherically symmetric ("s-type")
primitive Gaussians of the same orbital exponent, each centered on a
different atom. Closed forms for these integrals (in atomic units) are
standard (Szabo & Ostlund, *Modern Quantum Chemistry*, 1st ed. rev.,
Appendix A, eqs. A.9, A.11, A.33; Boys, *Proc. R. Soc. Lond. A* 200, 542
(1950), for the general Gaussian-product method) -- reproduced here in SI
units throughout, following chemistrykit's package-wide convention
(see :mod:`chemistrykit.constants`) of computing with actual physical
constants rather than adopting a rescaled (atomic) unit system.

Only s-type (l=0) Gaussians are implemented -- enough for a minimal 1s-like
basis on each hydrogen center -- so the general Boys function
:math:`F_n(x)` is only ever needed at :math:`n=0`.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from chemistrykit.constants import ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY

__all__ = ["GaussianPrimitive", "boys_f0", "overlap_integral", "kinetic_integral", "nuclear_attraction_integral"]

#: Coulomb's-law energy scale e^2/(4*pi*epsilon_0), in J*m -- the
#: recurring combination in every nuclear-attraction/electron-repulsion
#: integral below.
_COULOMB_CONSTANT = ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY)


class GaussianPrimitive:
    r"""A single normalized, s-type ("1s-like") Gaussian primitive, :math:`\chi(\mathbf r)=N e^{-\alpha|\mathbf r-\mathbf R|^2}`.

    The normalization constant :math:`N=(2\alpha/\pi)^{3/4}` is chosen so
    that :math:`\int|\chi|^2\,d^3r=1` (Szabo & Ostlund, Appendix A, eq.
    A.6, specialized to :math:`l=m=n=0`).

    Parameters
    ----------
    alpha : float
        Orbital exponent, in m^-2 (larger `alpha` -> a more sharply
        peaked, "tighter" orbital).
    center : array-like, shape (3,)
        Cartesian center of the Gaussian, in m.

    Examples
    --------
    >>> g = GaussianPrimitive(alpha=1.0e20, center=[0.0, 0.0, 0.0])
    >>> round(g.normalization, 6) == round((2.0 * 1.0e20 / np.pi) ** 0.75, 6)
    True
    """

    def __init__(self, alpha: float, center):
        if alpha <= 0:
            raise ValueError("alpha must be positive")
        self.alpha = float(alpha)
        self.center = np.asarray(center, dtype=np.float64)

    @property
    def normalization(self) -> float:
        """float: The normalization constant N."""
        return (2.0 * self.alpha / np.pi) ** 0.75


def boys_f0(x):
    r"""The zeroth-order Boys function, :math:`F_0(x)=\frac12\sqrt{\pi/x}\,\mathrm{erf}(\sqrt x)`.

    Arises from the Gaussian Coulomb (nuclear-attraction / electron-repulsion)
    integral (Boys, *Proc. R. Soc. Lond. A* 200, 542 (1950); Szabo &
    Ostlund, Appendix A, eq. A.32). Only :math:`n=0` is needed here since
    every basis function in this module is s-type.

    Parameters
    ----------
    x : float or array-like of float
        Must be non-negative.

    Returns
    -------
    float or ndarray
        :math:`F_0(0)=1` (the removable-singularity limit) is returned
        exactly at ``x=0``.

    Examples
    --------
    >>> round(float(boys_f0(0.0)), 6)
    1.0
    >>> import numpy as np
    >>> x = np.array([0.0, 1.0, 10.0])
    >>> bool(np.all(np.diff(boys_f0(x)) < 0))  # F0 is monotonically decreasing
    True
    """
    x = np.asarray(x, dtype=np.float64)
    scalar_input = x.ndim == 0
    x = np.atleast_1d(x)
    out = np.empty_like(x)
    small = x < 1.0e-12
    out[small] = 1.0
    xs = x[~small]
    out[~small] = 0.5 * np.sqrt(np.pi / xs) * erf(np.sqrt(xs))
    return float(out[0]) if scalar_input else out


def overlap_integral(a: GaussianPrimitive, b: GaussianPrimitive) -> float:
    r"""Overlap integral :math:`S_{ab}=\int\chi_a\chi_b\,d^3r` between two normalized s-Gaussians.

    .. math::

        S_{ab} = N_aN_b\left(\frac{\pi}{p}\right)^{3/2}
                 \exp\!\left(-\frac{\alpha_a\alpha_b}{p}|\mathbf R_a-\mathbf R_b|^2\right),
                 \qquad p=\alpha_a+\alpha_b

    (Szabo & Ostlund, Appendix A, eq. A.9, specialized to two s-functions
    of possibly different exponents.)

    Parameters
    ----------
    a, b : GaussianPrimitive

    Returns
    -------
    float
        Dimensionless (both `a` and `b` are normalized).

    Examples
    --------
    A Gaussian's overlap with itself is exactly 1 (it is normalized):

    >>> g = GaussianPrimitive(alpha=1.0e20, center=[0.0, 0.0, 0.0])
    >>> round(float(overlap_integral(g, g)), 9)
    1.0
    """
    p = a.alpha + b.alpha
    AB2 = float(np.sum((a.center - b.center) ** 2))
    return a.normalization * b.normalization * (np.pi / p) ** 1.5 * np.exp(-a.alpha * b.alpha / p * AB2)


def kinetic_integral(a: GaussianPrimitive, b: GaussianPrimitive) -> float:
    r"""Kinetic-energy integral :math:`T_{ab}=\int\chi_a\left(-\frac{\hbar^2}{2m_e}\nabla^2\right)\chi_b\,d^3r`.

    .. math::

        T_{ab} = N_aN_b\frac{\hbar^2}{m_e}\frac{\alpha_a\alpha_b}{p}
                 \left(3-\frac{2\alpha_a\alpha_b}{p}|\mathbf R_a-\mathbf R_b|^2\right)
                 \left(\frac{\pi}{p}\right)^{3/2}
                 \exp\!\left(-\frac{\alpha_a\alpha_b}{p}|\mathbf R_a-\mathbf R_b|^2\right)

    the standard atomic-units Gaussian kinetic-integral formula (Szabo &
    Ostlund, Appendix A, eq. A.11, there written for the operator
    :math:`-\frac12\nabla^2`, i.e. implicitly with :math:`\hbar=m_e=1`)
    restored to SI units by the explicit prefactor
    :math:`\hbar^2/m_e` in place of that implicit 1.

    Parameters
    ----------
    a, b : GaussianPrimitive

    Returns
    -------
    float
        Energy, in J.

    Examples
    --------
    The kinetic energy of a normalized Gaussian with itself is positive
    (a real physical kinetic energy):

    >>> g = GaussianPrimitive(alpha=1.0e20, center=[0.0, 0.0, 0.0])
    >>> bool(kinetic_integral(g, g) > 0)
    True
    """
    p = a.alpha + b.alpha
    AB2 = float(np.sum((a.center - b.center) ** 2))
    ab_over_p = a.alpha * b.alpha / p
    geometric = ab_over_p * (3.0 - 2.0 * ab_over_p * AB2) * (np.pi / p) ** 1.5 * np.exp(-ab_over_p * AB2)
    return a.normalization * b.normalization * (HBAR**2 / ELECTRON_MASS) * geometric


def nuclear_attraction_integral(a: GaussianPrimitive, b: GaussianPrimitive, Z: float, nucleus_center) -> float:
    r"""Nuclear-attraction integral :math:`V_{ab}=-Z\dfrac{e^2}{4\pi\varepsilon_0}\displaystyle\int\chi_a\frac{1}{|\mathbf r-\mathbf R_C|}\chi_b\,d^3r`.

    .. math::

        V_{ab} = -N_aN_b\,Z\frac{e^2}{4\pi\varepsilon_0}\frac{2\pi}{p}
                 \exp\!\left(-\frac{\alpha_a\alpha_b}{p}|\mathbf R_a-\mathbf R_b|^2\right)
                 F_0\!\left(p|\mathbf P-\mathbf R_C|^2\right)

    with :math:`\mathbf P=(\alpha_a\mathbf R_a+\alpha_b\mathbf R_b)/p` the
    Gaussian-product center and :math:`F_0` the Boys function
    (:func:`boys_f0`) (Szabo & Ostlund, Appendix A, eq. A.33, restored to
    SI units via the explicit Coulomb prefactor
    :math:`e^2/4\pi\varepsilon_0` in place of the atomic-units implicit 1).

    Parameters
    ----------
    a, b : GaussianPrimitive
    Z : float
        Nuclear charge, in units of the elementary charge (e.g. 1 for a proton).
    nucleus_center : array-like, shape (3,)
        Position of the attracting nucleus, in m.

    Returns
    -------
    float
        Energy, in J (negative -- attractive).

    Examples
    --------
    >>> g = GaussianPrimitive(alpha=1.0e20, center=[0.0, 0.0, 0.0])
    >>> bool(nuclear_attraction_integral(g, g, Z=1.0, nucleus_center=[0.0, 0.0, 0.0]) < 0)
    True
    """
    p = a.alpha + b.alpha
    AB2 = float(np.sum((a.center - b.center) ** 2))
    P = (a.alpha * a.center + b.alpha * b.center) / p
    nucleus_center = np.asarray(nucleus_center, dtype=np.float64)
    PC2 = float(np.sum((P - nucleus_center) ** 2))
    prefactor = a.normalization * b.normalization * Z * _COULOMB_CONSTANT * (2.0 * np.pi / p)
    return -prefactor * np.exp(-a.alpha * b.alpha / p * AB2) * boys_f0(p * PC2)
