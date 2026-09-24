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

from fractions import Fraction
from math import gcd, isinf, lcm

import numpy as np

__all__ = [
    "classify_crystal_system",
    "unit_cell_volume",
    "BRAVAIS_LATTICES",
    "cubic_lattice_points",
    "miller_indices_from_intercepts",
    "interplanar_angle_cubic",
]

#: dict: The 14 Bravais lattices (A. Bravais, 1848/1850), as the lattice
#: centerings allowed in each of the 7 crystal systems: ``"P"`` primitive,
#: ``"C"`` base-centered, ``"I"`` body-centered, ``"F"`` face-centered,
#: ``"R"`` rhombohedral (Ashcroft & Mermin, *Solid State Physics*, 1976,
#: Ch. 7). The counts sum to exactly 14.
BRAVAIS_LATTICES: dict = {
    "cubic": ("P", "I", "F"),
    "tetragonal": ("P", "I"),
    "orthorhombic": ("P", "C", "I", "F"),
    "hexagonal": ("P",),
    "trigonal": ("R",),
    "monoclinic": ("P", "C"),
    "triclinic": ("P",),
}

#: dict: Fractional lattice-point positions in the conventional cubic cell
#: for the three cubic Bravais centerings.
_CUBIC_CENTERINGS: dict = {
    "P": [(0.0, 0.0, 0.0)],
    "I": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    "F": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)],
}

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


def cubic_lattice_points(centering: str, a: float = 1.0, n_cells: int = 1) -> np.ndarray:
    r"""Cartesian lattice points of a cubic Bravais lattice inside an ``n_cells``-cube block of conventional cells.

    Generates every point :math:`a(\mathbf{n}+\mathbf{f})` with integer
    cell index :math:`\mathbf{n}` and centering offset :math:`\mathbf{f}`
    (see :data:`BRAVAIS_LATTICES`) that lies in the closed cube
    :math:`[0, n_{cells}a]^3` -- so corner and face points shared with
    neighboring cells are included, as in a textbook unit-cell drawing.

    Parameters
    ----------
    centering : {"P", "I", "F"}
        Primitive (simple), body-centered, or face-centered cubic.
    a : float, default 1.0
        Conventional cubic lattice constant.
    n_cells : int, default 1
        Number of conventional cells along each axis.

    Returns
    -------
    ndarray, shape (n_points, 3)

    Examples
    --------
    One conventional cell drawn with all shared corner/face points: 8
    corners (P), plus a body center (I), or plus 6 face centers (F):

    >>> [len(cubic_lattice_points(c)) for c in ("P", "I", "F")]
    [8, 9, 14]
    """
    offsets = _CUBIC_CENTERINGS[centering]
    eps = 1e-9
    points = []
    for i in range(n_cells + 1):
        for j in range(n_cells + 1):
            for k in range(n_cells + 1):
                for fx, fy, fz in offsets:
                    p = (i + fx, j + fy, k + fz)
                    if all(x <= n_cells + eps for x in p):
                        points.append(p)
    return a * np.array(points, dtype=np.float64)


def miller_indices_from_intercepts(a_intercept, b_intercept, c_intercept) -> tuple:
    r"""Miller indices :math:`(hkl)` of a crystal face from its axial intercepts (Hauy's law of rational indices).

    A face cutting the three crystallographic axes at :math:`p a`,
    :math:`q b`, :math:`r c` (intercepts in units of the cell edges) has
    indices proportional to :math:`(1/p, 1/q, 1/r)`, cleared of fractions
    and common factors to the smallest integer triple. Hauy's law states
    that for every natural face the intercept ratios are rational -- small
    whole-number ratios -- which is exactly what makes this integer triple
    exist (Hauy, *Traité de Minéralogie*, 1801; notation of W. H. Miller,
    *A Treatise on Crystallography*, 1839).

    Parameters
    ----------
    a_intercept, b_intercept, c_intercept : int, float, Fraction, or math.inf
        Intercepts in units of `a`, `b`, `c`; ``math.inf`` for a face
        parallel to that axis. Floats are converted with
        :meth:`fractions.Fraction.limit_denominator` (denominator <= 1000).
        Negative intercepts give negative indices.

    Returns
    -------
    tuple of (int, int, int)

    Examples
    --------
    A face cutting the axes at 1, 2, and 3 cell edges is the (632) face:

    >>> miller_indices_from_intercepts(1, 2, 3)
    (6, 3, 2)

    A cube face parallel to `b` and `c`:

    >>> import math
    >>> miller_indices_from_intercepts(1, math.inf, math.inf)
    (1, 0, 0)
    """
    reciprocals = []
    for x in (a_intercept, b_intercept, c_intercept):
        if isinstance(x, float) and isinf(x):
            reciprocals.append(Fraction(0))
            continue
        frac = x if isinstance(x, Fraction) else Fraction(x).limit_denominator(1000)
        if frac == 0:
            raise ValueError("a face cannot pass through the origin (zero intercept); shift the origin")
        reciprocals.append(1 / frac)
    if all(r == 0 for r in reciprocals):
        raise ValueError("at least one intercept must be finite")
    common_denominator = lcm(*(r.denominator for r in reciprocals))
    integers = [int(r * common_denominator) for r in reciprocals]
    divisor = gcd(*integers)
    return tuple(i // divisor for i in integers)


def interplanar_angle_cubic(hkl_1, hkl_2) -> float:
    r"""Angle between two crystal faces (planes) :math:`(h_1k_1l_1)` and :math:`(h_2k_2l_2)` of a cubic crystal, in degrees.

    .. math::

        \cos\phi = \frac{h_1h_2+k_1k_2+l_1l_2}
                         {\sqrt{h_1^2+k_1^2+l_1^2}\sqrt{h_2^2+k_2^2+l_2^2}}

    the angle between the face normals, which in a cubic crystal lie
    along :math:`[hkl]` (West, *Solid State Chemistry and its
    Applications*, 2nd ed., Ch. 1). It depends only on the indices, not
    on the lattice constant or the size of the faces -- Steno's constancy
    of interfacial angles.

    Parameters
    ----------
    hkl_1, hkl_2 : tuple of (int, int, int)
        Miller indices of the two planes (neither all zero).

    Returns
    -------
    float
        Angle between the plane normals, in degrees, in :math:`[0, 180]`.

    Examples
    --------
    Adjacent cube faces are perpendicular; a cube face and an octahedral
    face meet at :math:`\arccos(1/\sqrt3)\approx54.74°`:

    >>> interplanar_angle_cubic((1, 0, 0), (0, 1, 0))
    90.0
    >>> round(interplanar_angle_cubic((1, 0, 0), (1, 1, 1)), 2)
    54.74
    """
    v1 = np.asarray(hkl_1, dtype=np.float64)
    v2 = np.asarray(hkl_2, dtype=np.float64)
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        raise ValueError("(h, k, l) = (0, 0, 0) is not a plane")
    cos_phi = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_phi)))
