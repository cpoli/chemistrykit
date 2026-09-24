r"""Flory-Huggins lattice theory of polymer solutions.

P. J. Flory, *J. Chem. Phys.* 9, 660 (1941); M. L. Huggins, *J. Chem.
Phys.* 9, 440 (1941). See also Rubinstein & Colby, *Polymer Physics*
(2003), Ch. 4.

Polymer chains of :math:`N` segments and single-site solvent molecules
fill a lattice; with :math:`\phi` the polymer volume fraction, the free
energy of mixing per lattice site is

.. math::

    \frac{\Delta F_\text{mix}}{k_BT} = \frac{\phi}{N}\ln\phi
    + (1-\phi)\ln(1-\phi) + \chi\phi(1-\phi)

The first two (entropic) terms always favor mixing, but the polymer's
translational entropy is reduced by a factor :math:`1/N` because its
segments are connected; the dimensionless interaction parameter
:math:`\chi` captures the net energetic cost of polymer-solvent contacts.
The homogeneous solution becomes locally unstable (the spinodal) where
:math:`\partial^2\Delta F/\partial\phi^2=0`,

.. math::

    \chi_s(\phi) = \frac{1}{2}\left(\frac{1}{N\phi}+\frac{1}{1-\phi}\right)

whose minimum is the critical point
:math:`\phi_c=1/(1+\sqrt N)`, :math:`\chi_c=\tfrac12(1+1/\sqrt N)^2`.
As :math:`N\to\infty`, :math:`\chi_c\to1/2` -- the theta condition, where a
long chain's excluded volume :math:`v\propto1-2\chi` vanishes.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "flory_huggins_free_energy",
    "flory_huggins_spinodal_chi",
    "flory_huggins_critical_point",
]


def flory_huggins_free_energy(phi, N: float, chi: float):
    r"""Flory-Huggins free energy of mixing per lattice site, in units of :math:`k_BT`.

    Parameters
    ----------
    phi : float or array-like of float
        Polymer volume fraction, :math:`0<\phi<1`.
    N : float
        Degree of polymerization (lattice sites per chain).
    chi : float
        Flory-Huggins interaction parameter.

    Returns
    -------
    float or ndarray

    Examples
    --------
    With ``N=1`` and ``chi=0`` this is the ideal mixing entropy of a
    symmetric binary mixture, :math:`\ln\frac12` at :math:`\phi=1/2`:

    >>> round(float(flory_huggins_free_energy(0.5, N=1, chi=0.0)), 6)
    -0.693147
    """
    phi = np.asarray(phi, dtype=float)
    return phi / N * np.log(phi) + (1.0 - phi) * np.log(1.0 - phi) + chi * phi * (1.0 - phi)


def flory_huggins_spinodal_chi(phi, N: float):
    r"""Interaction parameter on the spinodal, :math:`\chi_s=\frac12[1/(N\phi)+1/(1-\phi)]`.

    Parameters
    ----------
    phi : float or array-like of float
        Polymer volume fraction.
    N : float
        Degree of polymerization.

    Returns
    -------
    float or ndarray

    Examples
    --------
    For a symmetric blend (``N=1``) the spinodal at :math:`\phi=1/2` is
    :math:`\chi=2`:

    >>> float(flory_huggins_spinodal_chi(0.5, N=1))
    2.0
    """
    phi = np.asarray(phi, dtype=float)
    return 0.5 * (1.0 / (N * phi) + 1.0 / (1.0 - phi))


def flory_huggins_critical_point(N: float) -> tuple[float, float]:
    r"""Critical point :math:`(\phi_c,\chi_c)` of a polymer solution.

    Parameters
    ----------
    N : float
        Degree of polymerization.

    Returns
    -------
    phi_c, chi_c : float
        :math:`\phi_c=1/(1+\sqrt N)` and :math:`\chi_c=\tfrac12(1+1/\sqrt N)^2`.

    Examples
    --------
    >>> phi_c, chi_c = flory_huggins_critical_point(100)
    >>> round(phi_c, 6), round(chi_c, 6)
    (0.090909, 0.605)

    For very long chains :math:`\chi_c` approaches the theta value 1/2:

    >>> round(flory_huggins_critical_point(1e12)[1], 5)
    0.5
    """
    s = np.sqrt(N)
    return float(1.0 / (1.0 + s)), float(0.5 * (1.0 + 1.0 / s) ** 2)
