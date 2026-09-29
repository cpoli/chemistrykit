r"""Minimal s-type Gaussian integrals, for the Hartree-Fock-style LCAO solvers.

Everything :mod:`chemistrykit.quantum.systems.hartree_fock` needs: the
overlap, kinetic-energy, one-electron nuclear-attraction and two-electron
repulsion integrals between normalized, spherically symmetric ("s-type")
primitive Gaussians centered on different atoms, and fixed contractions
of them (:class:`ContractedGaussian`, e.g. the STO-3G 1s function from
:func:`sto3g_1s`). Closed forms for these integrals (in atomic units) are
standard (Szabo & Ostlund, *Modern Quantum Chemistry*, 1st ed. rev.,
Appendix A, eqs. A.9, A.11, A.33, A.41; Boys, *Proc. R. Soc. Lond. A* 200,
542 (1950), for the general Gaussian-product method) -- reproduced here in
SI units throughout, following chemistrykit's package-wide convention
(see :mod:`chemistrykit.constants`) of computing with actual physical
constants rather than adopting a rescaled (atomic) unit system.

Only s-type (l=0) Gaussians are implemented -- enough for a minimal 1s
basis on each center of H2, HeH+ and similar first-row-free molecules --
so the general Boys function :math:`F_n(x)` is only ever needed at
:math:`n=0`.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from chemistrykit.constants import ELECTRON_MASS, ELEMENTARY_CHARGE, HBAR, VACUUM_PERMITTIVITY

__all__ = [
    "BOHR_RADIUS",
    "GaussianPrimitive",
    "ContractedGaussian",
    "sto3g_1s",
    "boys_f0",
    "overlap_integral",
    "kinetic_integral",
    "nuclear_attraction_integral",
    "electron_repulsion_integral",
    "contracted_one_electron",
    "contracted_electron_repulsion",
]

#: Coulomb's-law energy scale e^2/(4*pi*epsilon_0), in J*m -- the
#: recurring combination in every nuclear-attraction/electron-repulsion
#: integral below.
_COULOMB_CONSTANT = ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY)

#: The Bohr radius :math:`a_0=4\pi\varepsilon_0\hbar^2/(m_ee^2)`, in m -- the
#: length unit tabulated Gaussian exponents (in bohr^-2) are converted from.
BOHR_RADIUS = 4.0 * np.pi * VACUUM_PERMITTIVITY * HBAR**2 / (ELECTRON_MASS * ELEMENTARY_CHARGE**2)


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


def electron_repulsion_integral(a: GaussianPrimitive, b: GaussianPrimitive, c: GaussianPrimitive, d: GaussianPrimitive) -> float:
    r"""Two-electron repulsion integral :math:`(ab|cd)` between four normalized s-Gaussians.

    :math:`(ab|cd)=\frac{e^2}{4\pi\varepsilon_0}\iint\chi_a(1)\chi_b(1)\,r_{12}^{-1}\,\chi_c(2)\chi_d(2)\,d^3r_1d^3r_2`:

    .. math::

        (ab|cd) = N_aN_bN_cN_d\,\frac{e^2}{4\pi\varepsilon_0}
            \frac{2\pi^{5/2}}{pq\sqrt{p+q}}
            e^{-\frac{\alpha_a\alpha_b}{p}|\mathbf R_a-\mathbf R_b|^2
               -\frac{\alpha_c\alpha_d}{q}|\mathbf R_c-\mathbf R_d|^2}
            F_0\!\left(\frac{pq}{p+q}|\mathbf P-\mathbf Q|^2\right)

    with :math:`p=\alpha_a+\alpha_b`, :math:`q=\alpha_c+\alpha_d` and
    :math:`\mathbf P`, :math:`\mathbf Q` the two Gaussian-product centers
    (Szabo & Ostlund, Appendix A, eq. A.41, in chemists' notation, restored
    to SI units).

    Parameters
    ----------
    a, b, c, d : GaussianPrimitive
        `a`, `b` hold electron 1; `c`, `d` electron 2.

    Returns
    -------
    float
        Energy, in J (positive -- repulsive).

    Examples
    --------
    Two electrons in the same normalized Gaussian repel with
    :math:`(aa|aa)=\frac{e^2}{4\pi\varepsilon_0}\,2\sqrt{\alpha/\pi}`:

    >>> import numpy as np
    >>> g = GaussianPrimitive(alpha=1.0e20, center=[0.0, 0.0, 0.0])
    >>> expected = ELEMENTARY_CHARGE**2 / (4 * np.pi * VACUUM_PERMITTIVITY) * 2 * np.sqrt(1.0e20 / np.pi)
    >>> bool(np.isclose(electron_repulsion_integral(g, g, g, g), expected))
    True
    """
    p = a.alpha + b.alpha
    q = c.alpha + d.alpha
    AB2 = float(np.sum((a.center - b.center) ** 2))
    CD2 = float(np.sum((c.center - d.center) ** 2))
    P = (a.alpha * a.center + b.alpha * b.center) / p
    Q = (c.alpha * c.center + d.alpha * d.center) / q
    PQ2 = float(np.sum((P - Q) ** 2))
    norm = a.normalization * b.normalization * c.normalization * d.normalization
    geometric = 2.0 * np.pi**2.5 / (p * q * np.sqrt(p + q)) * np.exp(-a.alpha * b.alpha / p * AB2 - c.alpha * d.alpha / q * CD2)
    return norm * _COULOMB_CONSTANT * geometric * boys_f0(p * q / (p + q) * PQ2)


class ContractedGaussian:
    r"""A fixed linear combination of normalized s-type primitives on one center, :math:`\phi=\sum_k d_k\chi_k`.

    Parameters
    ----------
    exponents : array-like of float
        Primitive exponents, in m^-2.
    coefficients : array-like of float
        Contraction coefficients :math:`d_k` (for normalized primitives).
    center : array-like, shape (3,)
        Center, in m.

    Examples
    --------
    >>> phi = sto3g_1s(zeta=1.24, center=[0.0, 0.0, 0.0])
    >>> len(phi.primitives)
    3
    """

    def __init__(self, exponents, coefficients, center):
        exponents = np.asarray(exponents, dtype=np.float64)
        coefficients = np.asarray(coefficients, dtype=np.float64)
        if exponents.shape != coefficients.shape:
            raise ValueError("exponents and coefficients must have the same shape")
        self.center = np.asarray(center, dtype=np.float64)
        self.coefficients = coefficients
        self.primitives = [GaussianPrimitive(alpha, self.center) for alpha in exponents]


#: STO-3G least-squares fit of a zeta=1 Slater 1s orbital by three
#: Gaussians: exponents (bohr^-2) and contraction coefficients (W. J.
#: Hehre, R. F. Stewart & J. A. Pople, *J. Chem. Phys.* 51, 2657 (1969);
#: Szabo & Ostlund, eq. 3.225).
_STO3G_1S_EXPONENTS = np.array([0.109818, 0.405771, 2.22766])
_STO3G_1S_COEFFICIENTS = np.array([0.444635, 0.535328, 0.154329])


def sto3g_1s(zeta: float, center) -> ContractedGaussian:
    r"""The STO-3G contracted 1s function for Slater exponent `zeta` (exponents scale as :math:`\zeta^2`).

    Parameters
    ----------
    zeta : float
        Slater orbital exponent, in bohr^-1 (1.24 for H in molecules, 2.0925
        for He in HeH+; Szabo & Ostlund, Ch. 3.5.2).
    center : array-like, shape (3,)
        Center, in m.

    Returns
    -------
    ContractedGaussian

    Examples
    --------
    The contraction is normalized to within the fit's rounding:

    >>> phi = sto3g_1s(zeta=1.24, center=[0.0, 0.0, 0.0])
    >>> round(contracted_one_electron(overlap_integral, phi, phi), 5)
    1.0
    """
    return ContractedGaussian(_STO3G_1S_EXPONENTS * zeta**2 / BOHR_RADIUS**2, _STO3G_1S_COEFFICIENTS, center)


def contracted_one_electron(integral, a: ContractedGaussian, b: ContractedGaussian, *args) -> float:
    """Contract a primitive one-electron integral over two :class:`ContractedGaussian` functions.

    Parameters
    ----------
    integral : callable
        :func:`overlap_integral`, :func:`kinetic_integral` or
        :func:`nuclear_attraction_integral`.
    a, b : ContractedGaussian
    *args
        Extra arguments for `integral` (e.g. ``Z, nucleus_center``).

    Returns
    -------
    float
    """
    total = 0.0
    for da, pa in zip(a.coefficients, a.primitives, strict=True):
        for db, pb in zip(b.coefficients, b.primitives, strict=True):
            total += da * db * integral(pa, pb, *args)
    return float(total)


def contracted_electron_repulsion(a: ContractedGaussian, b: ContractedGaussian, c: ContractedGaussian, d: ContractedGaussian) -> float:
    """Contract :func:`electron_repulsion_integral` over four :class:`ContractedGaussian` functions.

    Parameters
    ----------
    a, b, c, d : ContractedGaussian

    Returns
    -------
    float
        :math:`(ab|cd)`, in J.
    """
    total = 0.0
    for da, pa in zip(a.coefficients, a.primitives, strict=True):
        for db, pb in zip(b.coefficients, b.primitives, strict=True):
            for dc, pc in zip(c.coefficients, c.primitives, strict=True):
                for dd, pd in zip(d.coefficients, d.primitives, strict=True):
                    total += da * db * dc * dd * electron_repulsion_integral(pa, pb, pc, pd)
    return float(total)
