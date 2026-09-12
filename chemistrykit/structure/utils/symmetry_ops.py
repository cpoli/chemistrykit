r"""Shared numerical machinery for testing candidate symmetry operations against a set of 3D coordinates.

Used by :mod:`chemistrykit.structure.systems.point_group` to find a
molecule's actual rotation axes, mirror planes, and inversion center by
direct geometric testing (see that module's docstring for the overall
approach) -- collected here, rather than inlined there, in the same
spirit as :mod:`chemistrykit.quantum.utils.secular_equation` hosting the
shared eigenvalue-solve machinery that
:mod:`chemistrykit.quantum.systems.huckel` and
:mod:`chemistrykit.quantum.systems.hartree_fock` both build a matrix for
and hand off to.
"""

from __future__ import annotations

import itertools

import numpy as np

__all__ = [
    "rotation_matrix",
    "reflection_matrix",
    "geometry_matches",
    "candidate_axes",
    "highest_rotation_order",
]


def rotation_matrix(axis, angle: float) -> np.ndarray:
    r"""Rodrigues' rotation-matrix formula for a right-handed rotation by `angle` about `axis`.

    .. math::

        R = I + \sin\theta\,[\hat n]_\times + (1-\cos\theta)[\hat n]_\times^2

    (O. Rodrigues, *J. Math. Pures Appl.* 5, 380 (1840); any standard
    classical-mechanics text, e.g. Goldstein, *Classical Mechanics*, 3rd
    ed., Ch. 4.7.)

    Parameters
    ----------
    axis : array-like of float, shape (3,)
        Rotation axis (need not be normalized).
    angle : float
        Rotation angle, in radians.

    Returns
    -------
    ndarray, shape (3, 3)

    Examples
    --------
    A 90-degree rotation about `z` sends `x` to `y`:

    >>> import numpy as np
    >>> R = rotation_matrix([0, 0, 1], np.pi / 2)
    >>> np.allclose(R @ [1, 0, 0], [0, 1, 0], atol=1e-10)
    True
    """
    axis = np.asarray(axis, dtype=np.float64)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1.0 - c
    return np.array(
        [
            [x * x * C + c, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, y * y * C + c, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, z * z * C + c],
        ]
    )


def reflection_matrix(normal) -> np.ndarray:
    r"""Householder reflection matrix for the plane through the origin with the given `normal`.

    .. math::

        R = I - 2\hat n\hat n^T

    Parameters
    ----------
    normal : array-like of float, shape (3,)
        Plane normal (need not be normalized).

    Returns
    -------
    ndarray, shape (3, 3)

    Examples
    --------
    Reflection through the `xy`-plane (normal along `z`) flips the sign
    of `z` only:

    >>> import numpy as np
    >>> R = reflection_matrix([0, 0, 1])
    >>> np.allclose(R @ [1, 2, 3], [1, 2, -3])
    True
    """
    normal = np.asarray(normal, dtype=np.float64)
    n = normal / np.linalg.norm(normal)
    return np.eye(3) - 2.0 * np.outer(n, n)


def geometry_matches(symbols, coords: np.ndarray, transformed: np.ndarray, tol: float) -> bool:
    """Whether `transformed` is a same-element permutation of `coords`, within distance `tol`.

    Parameters
    ----------
    symbols : sequence of str
        Element symbol of each atom, matched 1:1 with rows of `coords`.
    coords : ndarray, shape (n_atoms, 3)
        Reference (untransformed) coordinates.
    transformed : ndarray, shape (n_atoms, 3)
        Coordinates after applying a candidate symmetry operation.
    tol : float
        Maximum distance between a transformed atom and its matched
        original atom for the operation to count as a symmetry.

    Returns
    -------
    bool

    Examples
    --------
    >>> import numpy as np
    >>> coords = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    >>> geometry_matches(["H", "H"], coords, -coords, tol=1e-6)
    True
    >>> geometry_matches(["H", "F"], coords, -coords, tol=1e-6)
    False
    """
    n = len(symbols)
    used = [False] * n
    for i in range(n):
        found = False
        for j in range(n):
            if used[j] or symbols[j] != symbols[i]:
                continue
            if np.linalg.norm(transformed[i] - coords[j]) < tol:
                used[j] = True
                found = True
                break
        if not found:
            return False
    return True


def candidate_axes(coords: np.ndarray, tol: float = 1e-6) -> list:
    """Generate candidate symmetry-axis/mirror-normal directions from a set of atomic positions.

    Every non-central atom's position vector, every pairwise sum,
    difference, and cross product of two atoms' position vectors, and
    the inertia-tensor eigenvectors -- a broad, cheap candidate set that
    contains the true symmetry axis/mirror-plane normal for essentially
    any small, hand-built molecular geometry (see
    :mod:`chemistrykit.structure.systems.point_group`'s module docstring
    for the reasoning and its validated test cases).

    Parameters
    ----------
    coords : ndarray, shape (n_atoms, 3)
        Coordinates, already centered (e.g. at the centroid).
    tol : float, default 1e-6
        Vectors shorter than this are dropped (avoids a zero-vector
        candidate from an atom sitting exactly at the center).

    Returns
    -------
    list of ndarray, shape (3,)
        Unit vectors, deduplicated by direction (a vector and its
        negation are treated as the same candidate, since both a
        rotation axis and a mirror-plane normal are direction-agnostic).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n = len(coords)
    candidates = []
    for v in coords:
        if np.linalg.norm(v) > tol:
            candidates.append(v / np.linalg.norm(v))
    for i, j in itertools.combinations(range(n), 2):
        for combo in (coords[i] + coords[j], coords[i] - coords[j], np.cross(coords[i], coords[j])):
            if np.linalg.norm(combo) > tol:
                candidates.append(combo / np.linalg.norm(combo))
    inertia = np.zeros((3, 3))
    for v in coords:
        inertia += np.eye(3) * np.dot(v, v) - np.outer(v, v)
    _, eigenvectors = np.linalg.eigh(inertia)
    for k in range(3):
        candidates.append(eigenvectors[:, k])

    unique: list = []
    for c in candidates:
        if not any(abs(abs(np.dot(c, u)) - 1.0) < 1e-6 for u in unique):
            unique.append(c)
    return unique


def highest_rotation_order(symbols, coords: np.ndarray, axis, tol: float, max_n: int = 6) -> int:
    """Highest `n` (2..max_n) for which rotating `coords` by :math:`2\\pi/n` about `axis` is a symmetry operation.

    Parameters
    ----------
    symbols : sequence of str
    coords : ndarray, shape (n_atoms, 3)
        Already centered.
    axis : array-like of float, shape (3,)
        Candidate rotation axis.
    tol : float
        Forwarded to :func:`geometry_matches`.
    max_n : int, default 6
        Highest rotation order to test (sufficient for every point group
        except the linear ones, which
        :func:`chemistrykit.structure.systems.point_group.determine_point_group`
        special-cases before reaching this function).

    Returns
    -------
    int
        `1` if no rotation of order 2 through `max_n` about `axis` is a
        symmetry operation.
    """
    best = 1
    for n in range(2, max_n + 1):
        transformed = coords @ rotation_matrix(axis, 2.0 * np.pi / n).T
        if geometry_matches(symbols, coords, transformed, tol):
            best = n
    return best
