r"""The second virial coefficient from a pair potential (Mayer's cluster expansion).

J. E. Mayer, "The Statistical Mechanics of Condensing Systems. I,"
*J. Chem. Phys.* 5, 67-73 (1937), expanded an imperfect gas's
configurational integral in "clusters" of the Mayer function
:math:`f(r)=e^{-u(r)/k_BT}-1`, giving the virial expansion
:math:`PV_m/RT=1+B_2(T)/V_m+\cdots` with every coefficient an explicit
integral over the intermolecular potential. The first correction is

.. math::

    B_2(T) = -2\pi N_A\int_0^\infty\left(e^{-u(r)/k_BT}-1\right)r^2\,dr

(McQuarrie, *Statistical Mechanics*, Ch. 12). For hard spheres of
diameter :math:`\sigma` it is the temperature-independent
:math:`\frac{2\pi}{3}N_A\sigma^3`.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from scipy.integrate import quad

from chemistrykit.constants import K_B, NA

__all__ = ["second_virial_coefficient"]


def second_virial_coefficient(pair_potential: Callable[[float], float], T: float, length_scale: float, breakpoints: Sequence[float] = ()) -> float:
    r"""Return the molar second virial coefficient :math:`B_2(T)` of a spherical pair potential.

    Parameters
    ----------
    pair_potential : callable
        ``u(r)`` in J for a separation `r` in m; may return ``numpy.inf``
        inside a hard core.
    T : float
        Absolute temperature, in K.
    length_scale : float
        A characteristic molecular size (e.g. :math:`\sigma`), in m, used to
        split the integration range.
    breakpoints : sequence of float, optional
        Separations (m) where `u` is discontinuous (e.g. a square well's
        edges), passed on to the quadrature.

    Returns
    -------
    float
        :math:`B_2`, in m^3/mol (negative when attraction dominates).

    Examples
    --------
    Hard spheres give exactly :math:`\frac{2\pi}{3}N_A\sigma^3`:

    >>> sigma = 3.4e-10
    >>> hard_sphere = lambda r: np.inf if r < sigma else 0.0
    >>> B2 = second_virial_coefficient(hard_sphere, 300.0, sigma, breakpoints=[sigma])
    >>> round(B2 / (2 * np.pi / 3 * 6.02214076e23 * sigma**3), 8)
    1.0
    """
    if T <= 0 or length_scale <= 0:
        raise ValueError("T and length_scale must be positive")
    beta = 1.0 / (K_B * T)

    # Integrate in the reduced variable x = r / length_scale so the quadrature sees O(1) numbers.
    def integrand(x: float) -> float:
        if x <= 0.0:
            return 0.0
        with np.errstate(over="ignore", invalid="ignore"):
            u = float(pair_potential(x * length_scale))
        mayer_f = -1.0 if u == np.inf else np.expm1(-beta * u)
        return mayer_f * x * x

    edges = sorted({0.0, 0.5, 1.0, 2.0, 5.0, *(float(b) / length_scale for b in breakpoints if b > 0)})
    total = 0.0
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        total += quad(integrand, a, b, limit=500)[0]
    total += quad(integrand, edges[-1], np.inf, limit=500)[0]
    return -2.0 * np.pi * NA * length_scale**3 * total
