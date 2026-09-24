r"""Baeyer's angle-strain theory of carbocyclic rings.

Adolf von Baeyer argued in 1885 that a carbon atom's four bonds prefer the
tetrahedral angle of van't Hoff and Le Bel, :math:`\theta_T =
\arccos(-1/3) \approx 109.47^\circ`, and that a ring compound whose ring
is *planar* must bend its C-C-C angles to the polygon's interior angle,
:math:`180^\circ(n-2)/n`. Baeyer shared the distortion between the two
bonds at each carbon, so his angle strain per bond is

.. math::

    \delta(n) = \tfrac{1}{2}\left[\theta_T - \frac{180^\circ\,(n-2)}{n}\right]

(A. Baeyer, *Ber. Dtsch. Chem. Ges.* 18, 2269 (1885); see Eliel &
Wilen, *Stereochemistry of Organic Compounds* (Wiley, 1994), Ch. 11).
It predicts cyclopentane to be nearly strain-free and cyclopropane to be
badly strained, both correct. It also predicts growing strain for rings
of six or more, which is wrong: those rings pucker out of the plane
(Sachse, 1890) and can keep tetrahedral angles, as :func:`chair_cyclohexane_coordinates`
shows for the chair form of cyclohexane.
"""

from __future__ import annotations

import numpy as np

__all__ = ["TETRAHEDRAL_ANGLE", "planar_ring_angle", "baeyer_angle_strain", "chair_cyclohexane_coordinates"]

#: The ideal tetrahedral angle :math:`\arccos(-1/3)`, in degrees.
TETRAHEDRAL_ANGLE = float(np.degrees(np.arccos(-1.0 / 3.0)))


def planar_ring_angle(ring_size: int) -> float:
    """Interior angle of a regular planar polygon with `ring_size` vertices, in degrees.

    Parameters
    ----------
    ring_size : int
        Number of ring atoms, at least 3.

    Returns
    -------
    float

    Examples
    --------
    >>> planar_ring_angle(3), planar_ring_angle(6)
    (60.0, 120.0)
    """
    if ring_size < 3:
        raise ValueError("ring_size must be at least 3")
    return 180.0 * (ring_size - 2) / ring_size


def baeyer_angle_strain(ring_size: int) -> float:
    r"""Baeyer's angle strain per bond for a planar ring, in degrees.

    .. math::

        \delta(n) = \tfrac{1}{2}\left[109.47^\circ - \frac{180^\circ\,(n-2)}{n}\right]

    Parameters
    ----------
    ring_size : int
        Number of ring atoms, at least 3.

    Returns
    -------
    float
        Positive when the ring angle is squeezed below tetrahedral,
        negative when a planar ring would force it wider.

    Examples
    --------
    Baeyer's own values: 24 deg 44 min for cyclopropane, 9 deg 44 min for
    cyclobutane and only 0 deg 44 min for cyclopentane:

    >>> [round(baeyer_angle_strain(n), 2) for n in (3, 4, 5, 6)]
    [24.74, 9.74, 0.74, -5.26]
    """
    return 0.5 * (TETRAHEDRAL_ANGLE - planar_ring_angle(ring_size))


def chair_cyclohexane_coordinates(bond_length: float = 1.54) -> np.ndarray:
    r"""Carbon coordinates of an ideal chair cyclohexane with exactly tetrahedral C-C-C angles.

    The six carbons sit alternately at heights :math:`\pm z` on a circle
    of radius :math:`\rho`, 60 degrees apart. Requiring bond length `d`
    and bond angle :math:`\theta_T` fixes :math:`\rho^2 = 2d^2(1-\cos\theta_T)/3`
    (the 1-3 distance is :math:`\rho\sqrt3`) and then
    :math:`4z^2 = d^2 - \rho^2`.

    Parameters
    ----------
    bond_length : float, default 1.54
        C-C bond length, in angstrom.

    Returns
    -------
    ndarray, shape (6, 3)
        Ring atoms in order around the ring.

    Examples
    --------
    Every C-C-C angle of the puckered chair is tetrahedral, so Baeyer's
    predicted strain for a six-membered ring vanishes once the ring is not
    forced to be planar:

    >>> from chemistrykit.structure.core.base_system import angle_between
    >>> x = chair_cyclohexane_coordinates()
    >>> round(angle_between(x[0] - x[1], x[2] - x[1]), 4)
    109.4712
    """
    d = float(bond_length)
    rho = d * np.sqrt(2.0 * (1.0 - np.cos(np.radians(TETRAHEDRAL_ANGLE))) / 3.0)
    z = 0.5 * np.sqrt(d**2 - rho**2)
    phi = np.radians(60.0 * np.arange(6))
    heights = z * np.array([1, -1, 1, -1, 1, -1], dtype=np.float64)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), heights])
