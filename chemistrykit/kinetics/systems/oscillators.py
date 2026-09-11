"""The Brusselator: a minimal oscillating (limit-cycle) reaction network.

The kinetics analogue of the nonlinear oscillators studied elsewhere in
the *kit family (e.g. physicskit.chaos's continuous dynamical systems):
a chemical mechanism whose mass-action kinetics support a stable limit
cycle via a Hopf bifurcation, rather than settling to a fixed
concentration. Introduced by I. Prigogine & R. Lefever, *J. Chem. Phys.*
48, 1695 (1968); see also Nicolis & Prigogine, *Self-Organization in
Nonequilibrium Systems* (1977), Ch. 7.
"""

from __future__ import annotations

import numpy as np
from numba import njit
from numpy.typing import NDArray

from chemistrykit.kinetics.core.base_system import ReactionNetwork

__all__ = ["Brusselator"]


@njit(cache=True)
def _brusselator_rhs(state: NDArray[np.float64], t: float, params: NDArray[np.float64]) -> NDArray[np.float64]:
    """Brusselator vector field ``dX/dt, dY/dt = f(X, Y; A, B)``.

    Parameters
    ----------
    state : ndarray of float, shape (2,)
        State vector ``(X, Y)``.
    t : float
        Current time (unused; the system is autonomous).
    params : ndarray of float, shape (2,)
        Parameters ``(A, B)``.

    Returns
    -------
    ndarray of float, shape (2,)
    """
    A, B = params[0], params[1]
    X, Y = state[0], state[1]
    out = np.empty(2)
    out[0] = A - (B + 1.0) * X + X * X * Y
    out[1] = B * X - X * X * Y
    return out


class Brusselator(ReactionNetwork):
    r"""The Brusselator oscillating reaction network.

    Derived from four elementary steps (with B, D, E, and the source of
    X held as constant reservoirs; Prigogine & Lefever 1968):

    .. math::

        A &\to X \\
        2X + Y &\to 3X \\
        B + X &\to Y + D \\
        X &\to E

    which, in units where all four rate constants are 1, gives the net
    kinetics

    .. math::

        \frac{dX}{dt} = A - (B+1)X + X^2 Y, \qquad
        \frac{dY}{dt} = BX - X^2 Y

    The single fixed point :math:`(X^*, Y^*) = (A, B/A)` (see
    :meth:`fixed_point`) undergoes a Hopf bifurcation at :math:`B = 1 +
    A^2` (see :meth:`is_above_hopf_threshold`): below threshold the fixed
    point is a stable focus and concentrations relax to it; above
    threshold it is an unstable focus surrounded by a stable limit cycle,
    so concentrations settle into sustained oscillation instead.

    Parameters
    ----------
    X0, Y0 : float
        Initial concentrations.
    A, B : float
        Brusselator parameters (both must be positive).
    """

    species = ("X", "Y")

    def __init__(self, X0: float = 1.0, Y0: float = 1.0, A: float = 1.0, B: float = 3.0):
        if A <= 0 or B <= 0:
            raise ValueError("A and B must be positive")
        self.A = float(A)
        self.B = float(B)
        self.params = np.array([self.A, self.B])
        self._rhs_njit = _brusselator_rhs
        super().__init__([X0, Y0])

    def rhs(self, state, t):
        return np.asarray(_brusselator_rhs(np.asarray(state, dtype=np.float64), t, self.params))

    def fixed_point(self) -> np.ndarray:
        """Return the network's single steady state ``(X*, Y*) = (A, B/A)``.

        Returns
        -------
        ndarray, shape (2,)

        Examples
        --------
        The fixed point is exactly a zero of the vector field:

        >>> import numpy as np
        >>> system = Brusselator(A=1.0, B=3.0)
        >>> fp = system.fixed_point()
        >>> np.allclose(system.rhs(fp, 0.0), 0.0)
        True
        """
        return np.array([self.A, self.B / self.A])

    def is_above_hopf_threshold(self) -> bool:
        r"""Whether ``B > 1 + A^2``, the Hopf-bifurcation onset of sustained oscillation.

        Returns
        -------
        bool

        Examples
        --------
        >>> Brusselator(A=1.0, B=3.0).is_above_hopf_threshold()
        True
        >>> Brusselator(A=1.0, B=1.5).is_above_hopf_threshold()
        False
        """
        return self.B > 1.0 + self.A**2
