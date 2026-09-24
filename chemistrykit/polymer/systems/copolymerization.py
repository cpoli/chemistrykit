r"""Terminal-model copolymerization: the Mayo-Lewis copolymer equation.

F. R. Mayo & F. M. Lewis, *J. Am. Chem. Soc.* 66, 1594 (1944); see Odian,
*Principles of Polymerization*, 4th ed., Ch. 6.

If a growing chain's reactivity depends only on its terminal monomer
unit, four propagation steps with rate constants :math:`k_{11}, k_{12},
k_{21}, k_{22}` and a steady state in the two radical types give the
instantaneous mole fraction :math:`F_1` of monomer 1 entering the copolymer
from a feed of mole fraction :math:`f_1`:

.. math::

    F_1 = \frac{r_1f_1^2+f_1f_2}{r_1f_1^2+2f_1f_2+r_2f_2^2},
    \qquad r_1=\frac{k_{11}}{k_{12}},\; r_2=\frac{k_{22}}{k_{21}}

with :math:`f_2=1-f_1`. When both reactivity ratios are below (or both
above) one, the curve crosses the diagonal :math:`F_1=f_1` at the
azeotropic feed :math:`f_1^*=(1-r_2)/(2-r_1-r_2)`.
"""

from __future__ import annotations

import numpy as np

__all__ = ["mayo_lewis_copolymer_composition", "azeotropic_feed_composition"]


def mayo_lewis_copolymer_composition(f1, r1: float, r2: float):
    r"""Instantaneous copolymer composition :math:`F_1` from the Mayo-Lewis equation.

    Parameters
    ----------
    f1 : float or array-like of float
        Mole fraction of monomer 1 in the feed.
    r1, r2 : float
        Reactivity ratios.

    Returns
    -------
    float or ndarray

    Examples
    --------
    An ideal random copolymerization (:math:`r_1=r_2=1`) has
    :math:`F_1=f_1`:

    >>> float(mayo_lewis_copolymer_composition(0.3, 1.0, 1.0))
    0.3

    A perfectly alternating system (:math:`r_1=r_2=0`) gives
    :math:`F_1=1/2` at any feed:

    >>> float(mayo_lewis_copolymer_composition(0.9, 0.0, 0.0))
    0.5
    """
    f1 = np.asarray(f1, dtype=float)
    f2 = 1.0 - f1
    return (r1 * f1**2 + f1 * f2) / (r1 * f1**2 + 2.0 * f1 * f2 + r2 * f2**2)


def azeotropic_feed_composition(r1: float, r2: float) -> float:
    r"""Azeotropic feed :math:`f_1^*=(1-r_2)/(2-r_1-r_2)` at which :math:`F_1=f_1`.

    Parameters
    ----------
    r1, r2 : float
        Reactivity ratios (both < 1 or both > 1).

    Returns
    -------
    float

    Examples
    --------
    Styrene (1) / methyl methacrylate (2), :math:`r_1\approx0.52`,
    :math:`r_2\approx0.46`:

    >>> fstar = azeotropic_feed_composition(0.52, 0.46)
    >>> round(fstar, 4)
    0.5294
    >>> round(float(mayo_lewis_copolymer_composition(fstar, 0.52, 0.46)), 4)
    0.5294
    """
    return (1.0 - r2) / (2.0 - r1 - r2)
