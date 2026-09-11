r"""Periodic boundary conditions: minimum-image displacement and wrapping.

See Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch.
1.5.2 for the minimum-image convention itself.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = ["minimum_image_displacement", "wrap_positions"]


def minimum_image_displacement(r_i, r_j, box_length: Optional[float]):
    r"""Displacement :math:`\vec r_i - \vec r_j` under the minimum-image convention.

    .. math::

        \Delta\vec r = (\vec r_i - \vec r_j)
        - L \, \mathrm{round}\!\left(\frac{\vec r_i - \vec r_j}{L}\right)

    for a cubic periodic box of side `L` -- the displacement to the
    *closest* periodic image of particle `j` relative to particle `i`,
    valid whenever `L` is at least twice the interaction cutoff (Allen &
    Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch. 1.5.2). If
    `box_length` is ``None`` or ``<= 0``, plain (non-periodic)
    displacement is returned unchanged.

    Parameters
    ----------
    r_i, r_j : array-like
        Broadcastable position array(s).
    box_length : float or None

    Returns
    -------
    ndarray

    Examples
    --------
    Two particles near opposite edges of a box of length 10 are actually
    close together once the box wraps around:

    >>> import numpy as np
    >>> r_i = np.array([0.5, 0.0, 0.0])
    >>> r_j = np.array([9.5, 0.0, 0.0])
    >>> minimum_image_displacement(r_i, r_j, box_length=10.0)
    array([1., 0., 0.])

    With no box, this is just the ordinary displacement:

    >>> minimum_image_displacement(r_i, r_j, box_length=None)
    array([-9.,  0.,  0.])
    """
    r_i = np.asarray(r_i, dtype=np.float64)
    r_j = np.asarray(r_j, dtype=np.float64)
    diff = r_i - r_j
    if box_length is not None and box_length > 0:
        diff = diff - box_length * np.round(diff / box_length)
    return diff


def wrap_positions(positions, box_length: Optional[float]):
    """Wrap positions into the primary periodic cell ``[0, box_length)``.

    Parameters
    ----------
    positions : array-like
    box_length : float or None
        If ``None`` or ``<= 0``, `positions` is returned unchanged (as a
        float array).

    Returns
    -------
    ndarray

    Examples
    --------
    >>> import numpy as np
    >>> wrap_positions(np.array([-0.5, 10.2, 5.0]), box_length=10.0)
    array([9.5, 0.2, 5. ])
    """
    positions = np.asarray(positions, dtype=np.float64)
    if box_length is None or box_length <= 0:
        return positions
    return np.mod(positions, box_length)
