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

__all__ = ["Brusselator", "Oregonator"]


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


@njit(cache=True)
def _oregonator_rhs(state: NDArray[np.float64], t: float, params: NDArray[np.float64]) -> NDArray[np.float64]:
    """Scaled three-variable Oregonator vector field.

    Parameters
    ----------
    state : ndarray of float, shape (3,)
        Scaled state ``(x, y, z)`` (HBrO2, Br-, oxidized catalyst).
    t : float
        Current time (unused; the system is autonomous).
    params : ndarray of float, shape (4,)
        Parameters ``(epsilon, epsilon_prime, q, f)``.

    Returns
    -------
    ndarray of float, shape (3,)
    """
    eps, eps_p, q, f = params[0], params[1], params[2], params[3]
    x, y, z = state[0], state[1], state[2]
    out = np.empty(3)
    out[0] = (q * y - x * y + x * (1.0 - x)) / eps
    out[1] = (-q * y - x * y + f * z) / eps_p
    out[2] = x - z
    return out


class Oregonator(ReactionNetwork):
    r"""The Field-Noyes Oregonator model of the Belousov-Zhabotinsky reaction.

    Field and Noyes (*J. Chem. Phys.* 60, 1877 (1974)) reduced the FKN
    mechanism of the cerium-catalyzed bromate/malonic-acid (BZ) reaction
    to five steps among three intermediates -- X = HBrO2, Y = Br-, and
    Z = Ce(IV) -- with the bromate (A) and organic substrate (B) pools
    held constant:

    .. math::

        A + Y &\to X + P, \qquad X + Y \to 2P, \qquad A + X \to 2X + 2Z, \\
        2X &\to A + P, \qquad B + Z \to f\,Y

    In the standard dimensionless scaling (Tyson, in R. J. Field and M.
    Burger (eds.), *Oscillations and Traveling Waves in Chemical Systems*,
    Wiley 1985, Ch. 3) the kinetics read

    .. math::

        \epsilon \frac{dx}{dt} = qy - xy + x(1-x), \quad
        \epsilon' \frac{dy}{dt} = -qy - xy + fz, \quad
        \frac{dz}{dt} = x - z

    The system is stiff (:math:`\epsilon' \ll \epsilon \ll 1`), so
    integrate with ``method="dopri5"``. For the default parameters it
    relaxes onto a stable limit cycle (sustained BZ oscillation); for
    large stoichiometric factor `f` (e.g. ``f=3``) its steady state is
    stable instead.

    Parameters
    ----------
    x0, y0, z0 : float
        Initial scaled concentrations.
    epsilon, epsilon_prime, q, f : float
        Scaled Oregonator parameters (all positive).
    """

    species = ("X", "Y", "Z")

    def __init__(
        self,
        x0: float = 0.5,
        y0: float = 0.1,
        z0: float = 0.1,
        epsilon: float = 4e-2,
        epsilon_prime: float = 4e-4,
        q: float = 8e-4,
        f: float = 1.0,
    ):
        if min(epsilon, epsilon_prime, q, f) <= 0:
            raise ValueError("epsilon, epsilon_prime, q and f must be positive")
        self.epsilon = float(epsilon)
        self.epsilon_prime = float(epsilon_prime)
        self.q = float(q)
        self.f = float(f)
        self.params = np.array([self.epsilon, self.epsilon_prime, self.q, self.f])
        self._rhs_njit = _oregonator_rhs
        super().__init__([x0, y0, z0])

    def rhs(self, state, t):
        return np.asarray(_oregonator_rhs(np.asarray(state, dtype=np.float64), t, self.params))

    def fixed_point(self) -> np.ndarray:
        r"""Return the positive steady state ``(x*, y*, z*)``.

        Setting all three derivatives to zero gives :math:`z^* = x^*`,
        :math:`y^* = f x^*/(q + x^*)`, and :math:`x^*` as the positive
        root of :math:`x^2 - (1 - f - q)x - q(1 + f) = 0`.

        Returns
        -------
        ndarray, shape (3,)

        Examples
        --------
        >>> import numpy as np
        >>> system = Oregonator(f=1.0)
        >>> np.allclose(system.rhs(system.fixed_point(), 0.0), 0.0, atol=1e-9)
        True
        """
        b = 1.0 - self.f - self.q
        x = 0.5 * (b + np.sqrt(b * b + 4.0 * self.q * (1.0 + self.f)))
        return np.array([x, self.f * x / (self.q + x), x])
