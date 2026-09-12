r"""Ideal random-walk chain statistics, and the Flory exponent for real chains.

See Rubinstein & Colby, *Polymer Physics* (2003), Ch. 2 ("Ideal chains")
and Ch. 3 ("Real chains"), and P. J. Flory, *Principles of Polymer
Chemistry* (1953), Ch. X, for the original scaling argument.

An **ideal chain** models a polymer as a simple random walk of :math:`n`
freely-jointed segments (Kuhn segments) of length :math:`b`, ignoring
excluded volume (the chain can pass through itself) -- exactly realized
physically by a real chain at its *theta temperature*, where segment-
segment and segment-solvent interactions happen to cancel. Its two
headline results are exact consequences of the central-limit theorem
applied to a sum of :math:`n` random steps:

.. math::

    \langle R^2\rangle = nb^2, \qquad \langle R_g^2\rangle = \frac{nb^2}{6}

A **real chain** additionally accounts for excluded volume, and its size
scales instead as :math:`R \sim bn^\nu` with a solvent-quality-dependent
Flory exponent :math:`\nu` (Flory's mean-field estimate, later refined by
renormalization-group theory to :math:`\nu\approx0.588` for a good
solvent in 3D -- the value used here, per Rubinstein & Colby Table 3.2):

* **theta solvent** (:math:`\nu=1/2`): excluded-volume and solvent
  interactions cancel; the chain behaves exactly as an ideal chain.
* **good solvent** (:math:`\nu\approx3/5`): the solvent favors chain-
  solvent contacts, swelling the chain relative to ideal.
* **poor solvent** (:math:`\nu=1/3`): the solvent favors chain-chain
  contacts, collapsing the chain to a dense globule of essentially
  constant *density* (: math:`R^3\sim n`, hence :math:`\nu=1/3`).

**Approximation flagged explicitly**: for a *real* chain (any solvent
other than theta), the radius-of-gyration-to-end-to-end-distance ratio
used below (:math:`\langle R_g^2\rangle = \langle R^2\rangle/6`) is
carried over from the exact ideal-chain result as a scaling-law estimate
-- the true prefactor for a self-avoiding walk is solvent- and model-
dependent and not exactly 1/6 (Rubinstein & Colby, *op. cit.*, Ch. 3.2);
only the ideal-chain (theta-solvent) value is exact.
"""

from __future__ import annotations

from chemistrykit.polymer.core.base_system import PolymerChainModel

__all__ = ["FLORY_EXPONENTS", "flory_exponent", "IdealChain", "RealChain"]

#: dict: Standard Flory exponents nu for R ~ b*n**nu, by solvent quality
#: (Rubinstein & Colby, *Polymer Physics*, 2003, Table 3.2).
FLORY_EXPONENTS = {
    "theta": 0.5,
    "good": 3.0 / 5.0,
    "poor": 1.0 / 3.0,
}


def flory_exponent(solvent: str) -> float:
    r"""Return the standard Flory scaling exponent :math:`\nu` for a named solvent quality.

    Parameters
    ----------
    solvent : {"theta", "good", "poor"}
        Solvent quality (see the module docstring).

    Returns
    -------
    float

    Examples
    --------
    >>> flory_exponent("theta")
    0.5
    >>> round(flory_exponent("good"), 3)
    0.6
    >>> round(flory_exponent("poor"), 4)
    0.3333
    """
    try:
        return FLORY_EXPONENTS[solvent]
    except KeyError:
        raise ValueError(f"unknown solvent quality {solvent!r}; expected one of {sorted(FLORY_EXPONENTS)}") from None


class IdealChain(PolymerChainModel):
    r"""An ideal (theta-solvent) polymer chain: an exact random walk of `n` segments.

    Examples
    --------
    The mean-square end-to-end distance is *exactly* :math:`nb^2`
    (Rubinstein & Colby, *Polymer Physics*, 2003, eq. 2.44):

    >>> chain = IdealChain()
    >>> chain.mean_square_end_to_end(n=100, b=0.5)
    25.0

    The radius of gyration is exactly :math:`R/\sqrt{6}`:

    >>> import numpy as np
    >>> n, b = 200, 0.6
    >>> ratio = chain.mean_square_end_to_end(n, b) / chain.mean_square_radius_of_gyration(n, b)
    >>> round(float(ratio), 6)
    6.0

    Doubling `n` exactly doubles :math:`\langle R^2\rangle` (linear
    scaling, :math:`\nu=1/2` since :math:`R^2\sim n^{2\nu}=n`):

    >>> chain.mean_square_end_to_end(n=400, b=0.5) / chain.mean_square_end_to_end(n=200, b=0.5)
    2.0
    """

    def mean_square_end_to_end(self, n, b):
        return n * b**2

    def mean_square_radius_of_gyration(self, n, b):
        return n * b**2 / 6.0


class RealChain(PolymerChainModel):
    r"""A real (excluded-volume) polymer chain: :math:`R = bn^\nu` with Flory exponent `nu`.

    Parameters
    ----------
    nu : float
        Flory scaling exponent. Use :meth:`theta_solvent`,
        :meth:`good_solvent`, or :meth:`poor_solvent` for the standard
        named values rather than passing `nu` directly, unless modeling a
        nonstandard solvent quality.

    Examples
    --------
    A real chain in a good solvent is more swollen (larger) than an
    ideal chain of the same `n`:

    >>> ideal = IdealChain()
    >>> good = RealChain.good_solvent()
    >>> n, b = 1000, 0.5
    >>> bool(good.end_to_end_distance(n, b) > ideal.end_to_end_distance(n, b))
    True

    A real chain in a poor solvent is more collapsed (smaller):

    >>> poor = RealChain.poor_solvent()
    >>> bool(poor.end_to_end_distance(n, b) < ideal.end_to_end_distance(n, b))
    True

    At the theta point, :class:`RealChain` reduces exactly to the ideal
    chain's end-to-end distance (though not, by construction here, to its
    exact radius-of-gyration prefactor -- see the module docstring):

    >>> theta = RealChain.theta_solvent()
    >>> round(float(theta.end_to_end_distance(n, b)), 6) == round(float(ideal.end_to_end_distance(n, b)), 6)
    True
    """

    def __init__(self, nu: float):
        if not (0.0 < nu < 1.0):
            raise ValueError("nu must be in (0, 1)")
        self.nu = float(nu)

    @classmethod
    def theta_solvent(cls) -> RealChain:
        """Build a :class:`RealChain` with the theta-solvent Flory exponent (nu=1/2)."""
        return cls(FLORY_EXPONENTS["theta"])

    @classmethod
    def good_solvent(cls) -> RealChain:
        """Build a :class:`RealChain` with the good-solvent Flory exponent (nu=3/5)."""
        return cls(FLORY_EXPONENTS["good"])

    @classmethod
    def poor_solvent(cls) -> RealChain:
        """Build a :class:`RealChain` with the poor-solvent Flory exponent (nu=1/3)."""
        return cls(FLORY_EXPONENTS["poor"])

    def mean_square_end_to_end(self, n, b):
        return (b * n**self.nu) ** 2

    def mean_square_radius_of_gyration(self, n, b):
        """See the module docstring: the 1/6 prefactor is an ideal-chain approximation, not exact for real chains."""
        return self.mean_square_end_to_end(n, b) / 6.0
