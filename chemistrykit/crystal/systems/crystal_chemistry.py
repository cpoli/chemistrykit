r"""Crystal-chemistry rules of thumb from ionic radii: Pauling's radius-ratio rule and Goldschmidt's tolerance factor.

Both treat ions as hard spheres of fixed radius (e.g. Shannon's effective
ionic radii, :data:`chemistrykit.crystal.utils.reference_data.SHANNON_IONIC_RADII_PM`)
and ask purely geometric questions: how many anions fit around a cation
while still touching it (Pauling), and how well an :math:`ABX_3` cation
pair fits the ideal cubic perovskite framework (Goldschmidt). See West,
*Solid State Chemistry and its Applications*, 2nd ed. (2014), Ch. 3, for
both.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["RadiusRatioPrediction", "radius_ratio_coordination", "goldschmidt_tolerance_factor", "RADIUS_RATIO_LIMITS"]

#: tuple of (float, int, str): Lower radius-ratio limit, coordination
#: number, and coordination polyhedron, in increasing order. Each limit is
#: the exact ratio :math:`r_+/r_-` at which anions touching the cation
#: also touch each other: :math:`\sqrt{4/3}-1` (triangle),
#: :math:`\sqrt{3/2}-1` (tetrahedron), :math:`\sqrt2-1` (octahedron),
#: :math:`\sqrt3-1` (cube).
RADIUS_RATIO_LIMITS: tuple = (
    (0.0, 2, "linear"),
    (float(np.sqrt(4.0 / 3.0) - 1.0), 3, "trigonal planar"),
    (float(np.sqrt(1.5) - 1.0), 4, "tetrahedral"),
    (float(np.sqrt(2.0) - 1.0), 6, "octahedral"),
    (float(np.sqrt(3.0) - 1.0), 8, "cubic"),
)


@dataclass
class RadiusRatioPrediction:
    """The coordination predicted by the radius-ratio rule for one cation-anion pair."""

    radius_ratio: float
    """float: :math:`r_+/r_-`."""

    coordination_number: int
    """int: Predicted number of anions around the cation."""

    geometry: str
    """str: Predicted coordination polyhedron."""


def radius_ratio_coordination(r_cation: float, r_anion: float) -> RadiusRatioPrediction:
    r"""Predict a cation's coordination number from the radius ratio :math:`r_+/r_-` (Pauling's first rule).

    Pauling's first rule (L. Pauling, *J. Am. Chem. Soc.* 51 (1929),
    1010-1026): a coordination polyhedron of anions forms around each
    cation, with the coordination number fixed by the radius ratio. The
    largest polyhedron is chosen whose anions can all touch the cation
    without overlapping one another, i.e. the highest coordination whose
    lower limit in :data:`RADIUS_RATIO_LIMITS` does not exceed
    :math:`r_+/r_-`. For example, octahedral coordination needs
    :math:`r_+/r_-\ge\sqrt2-1\approx0.414`, since at that ratio four anions
    around the cation's equator touch each other.

    Parameters
    ----------
    r_cation, r_anion : float
        Ionic radii, in the same unit.

    Returns
    -------
    RadiusRatioPrediction

    Examples
    --------
    NaCl (Shannon radii 102 and 181 pm) is predicted octahedral, as in
    rock salt:

    >>> p = radius_ratio_coordination(102.0, 181.0)
    >>> p.coordination_number, p.geometry
    (6, 'octahedral')

    CsCl's large cation gives eightfold (cubic) coordination:

    >>> radius_ratio_coordination(174.0, 181.0).coordination_number
    8
    """
    if r_cation <= 0 or r_anion <= 0:
        raise ValueError("ionic radii must be positive")
    ratio = r_cation / r_anion
    chosen = RADIUS_RATIO_LIMITS[0]
    for limit in RADIUS_RATIO_LIMITS:
        if ratio >= limit[0]:
            chosen = limit
    return RadiusRatioPrediction(radius_ratio=float(ratio), coordination_number=chosen[1], geometry=chosen[2])


def goldschmidt_tolerance_factor(r_a, r_b, r_x):
    r"""Goldschmidt's tolerance factor for an :math:`ABX_3` perovskite.

    .. math::

        t = \frac{r_A + r_X}{\sqrt{2}\,(r_B + r_X)}

    In the ideal cubic perovskite the A-X distance is the half face
    diagonal, :math:`a/\sqrt2`, and the B-X distance is the half edge,
    :math:`a/2`, so touching hard spheres give exactly :math:`t=1`
    (V. M. Goldschmidt, *Naturwissenschaften* 14 (1926), 477-485).
    Roughly, :math:`0.9\lesssim t\lesssim1` gives a cubic perovskite,
    :math:`t` somewhat below 0.9 a tilted, lower-symmetry perovskite,
    and :math:`t>1` a hexagonal or tetragonal (ferroelectric) distortion.

    Parameters
    ----------
    r_a, r_b, r_x : float or array-like of float
        Ionic radii of the A cation (12-coordinate), B cation
        (6-coordinate), and X anion, in the same unit.

    Returns
    -------
    float or ndarray

    Examples
    --------
    SrTiO3, the archetypal cubic perovskite (Shannon radii Sr2+ XII
    144 pm, Ti4+ VI 60.5 pm, O2- 140 pm), is almost exactly ideal:

    >>> round(goldschmidt_tolerance_factor(144.0, 60.5, 140.0), 3)
    1.002

    Radii that satisfy :math:`r_A+r_X=\sqrt2(r_B+r_X)` give :math:`t=1` exactly:

    >>> import numpy as np
    >>> r_b, r_x = 60.0, 140.0
    >>> round(goldschmidt_tolerance_factor(np.sqrt(2.0) * (r_b + r_x) - r_x, r_b, r_x), 12)
    1.0
    """
    r_a = np.asarray(r_a, dtype=np.float64)
    r_b = np.asarray(r_b, dtype=np.float64)
    r_x = np.asarray(r_x, dtype=np.float64)
    t = (r_a + r_x) / (np.sqrt(2.0) * (r_b + r_x))
    return float(t) if t.ndim == 0 else t
