r"""The 7 crystal systems, and the general unit-cell-volume formula.

See Ashcroft & Mermin, *Solid State Physics* (1976), Ch. 7, or West,
*Solid State Chemistry and its Applications*, 2nd ed. (2014), Ch. 1, for
the standard classification of the 7 crystal systems by the symmetry
constraints they impose on the six lattice parameters
:math:`(a,b,c,\alpha,\beta,\gamma)` (edge lengths and interaxial angles
of the conventional unit cell):

======================  ===============================  ==============================
System                  Edge-length constraint            Angle constraint
======================  ===============================  ==============================
cubic                   :math:`a=b=c`                     :math:`\alpha=\beta=\gamma=90°`
tetragonal              :math:`a=b\neq c`                 :math:`\alpha=\beta=\gamma=90°`
orthorhombic            :math:`a\neq b\neq c`             :math:`\alpha=\beta=\gamma=90°`
hexagonal               :math:`a=b\neq c`                 :math:`\alpha=\beta=90°,\gamma=120°`
trigonal (rhombohedral) :math:`a=b=c`                     :math:`\alpha=\beta=\gamma\neq90°`
monoclinic              :math:`a\neq b\neq c`             :math:`\alpha=\gamma=90°\neq\beta`
triclinic               :math:`a\neq b\neq c`             :math:`\alpha\neq\beta\neq\gamma`
======================  ===============================  ==============================

Unlike :func:`chemistrykit.structure.systems.point_group.determine_point_group`,
which tests candidate symmetry *operations* against a molecule's actual 3D
coordinates, classifying a crystal system here is purely a lookup on
equalities/inequalities among the six lattice parameters themselves --
see the note in :mod:`chemistrykit.crystal.core.base_system` on why the
two modules do not share code despite the conceptual overlap.
"""

from __future__ import annotations

import numpy as np

__all__ = ["classify_crystal_system", "unit_cell_volume"]

_CRYSTAL_SYSTEMS = (
    "cubic",
    "tetragonal",
    "orthorhombic",
    "hexagonal",
    "trigonal",
    "monoclinic",
    "triclinic",
)


def classify_crystal_system(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, tol: float = 1e-6) -> str:
    r"""Classify a unit cell into one of the 7 crystal systems from its lattice parameters.

    Parameters
    ----------
    a, b, c : float
        Unit-cell edge lengths (any consistent length unit).
    alpha, beta, gamma : float
        Interaxial angles, in degrees (:math:`\alpha` between `b` and
        `c`, :math:`\beta` between `a` and `c`, :math:`\gamma` between
        `a` and `b`, the standard crystallographic convention).
    tol : float, default 1e-6
        Absolute tolerance for treating two lengths, or an angle and a
        reference angle, as equal.

    Returns
    -------
    str
        One of ``"cubic"``, ``"tetragonal"``, ``"orthorhombic"``,
        ``"hexagonal"``, ``"trigonal"``, ``"monoclinic"``, ``"triclinic"``.

    Examples
    --------
    A cell with all edges and angles equal to a right angle is cubic:

    >>> classify_crystal_system(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    'cubic'

    Hexagonal graphite-like cell (:math:`a=b\neq c`, :math:`\gamma=120°`):

    >>> classify_crystal_system(2.46, 2.46, 6.71, 90.0, 90.0, 120.0)
    'hexagonal'

    A cell with no special equalities at all is triclinic:

    >>> classify_crystal_system(5.0, 6.0, 7.0, 80.0, 85.0, 95.0)
    'triclinic'
    """

    def eq(x: float, y: float) -> bool:
        return abs(x - y) < tol

    ab, bc, ac = eq(a, b), eq(b, c), eq(a, c)
    all_edges_equal = ab and bc
    no_edges_equal = not ab and not bc and not ac

    right_angles = eq(alpha, 90.0) and eq(beta, 90.0) and eq(gamma, 90.0)

    if all_edges_equal and right_angles:
        return "cubic"
    if all_edges_equal and eq(alpha, beta) and eq(beta, gamma) and not eq(alpha, 90.0):
        return "trigonal"
    if ab and not bc and eq(alpha, 90.0) and eq(beta, 90.0) and eq(gamma, 120.0):
        return "hexagonal"
    if ab and not bc and right_angles:
        return "tetragonal"
    if no_edges_equal and right_angles:
        return "orthorhombic"
    if no_edges_equal and eq(alpha, 90.0) and eq(gamma, 90.0) and not eq(beta, 90.0):
        return "monoclinic"
    return "triclinic"


def unit_cell_volume(a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
    r"""General unit-cell volume from the six lattice parameters, valid for any crystal system.

    .. math::

        V = abc\sqrt{1-\cos^2\alpha-\cos^2\beta-\cos^2\gamma+2\cos\alpha\cos\beta\cos\gamma}

    (Ashcroft & Mermin, *Solid State Physics*, 1976, derived from the
    scalar triple product of the three edge vectors). Reduces to
    :math:`V=abc` for any right-angle cell (cubic/tetragonal/
    orthorhombic) and to :math:`V=\frac{\sqrt3}{2}a^2c` for hexagonal.

    Parameters
    ----------
    a, b, c : float or array-like of float
        Unit-cell edge lengths.
    alpha, beta, gamma : float or array-like of float
        Interaxial angles, in degrees.

    Returns
    -------
    float or ndarray

    Examples
    --------
    A right-angle (orthorhombic/cubic/tetragonal) cell reduces to :math:`abc`:

    >>> round(float(unit_cell_volume(2.0, 3.0, 4.0, 90.0, 90.0, 90.0)), 6)
    24.0

    A hexagonal cell (:math:`a=b`, :math:`\gamma=120°`) matches the
    textbook shortcut :math:`\frac{\sqrt3}{2}a^2c`:

    >>> import numpy as np
    >>> a, c = 2.46, 6.71
    >>> V_general = unit_cell_volume(a, a, c, 90.0, 90.0, 120.0)
    >>> V_hex_shortcut = (np.sqrt(3.0) / 2.0) * a**2 * c
    >>> round(float(V_general), 6) == round(float(V_hex_shortcut), 6)
    True
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64)
    ca = np.cos(np.radians(np.asarray(alpha, dtype=np.float64)))
    cb = np.cos(np.radians(np.asarray(beta, dtype=np.float64)))
    cg = np.cos(np.radians(np.asarray(gamma, dtype=np.float64)))
    factor = np.sqrt(np.clip(1.0 - ca**2 - cb**2 - cg**2 + 2.0 * ca * cb * cg, 0.0, None))
    result = a * b * c * factor
    return float(result) if result.ndim == 0 else result
