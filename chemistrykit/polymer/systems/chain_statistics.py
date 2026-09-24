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
Flory exponent :math:`\nu`:

* **theta solvent** (:math:`\nu=1/2`): excluded-volume and solvent
  interactions cancel; the chain behaves exactly as an ideal chain.
* **good solvent** (:math:`\nu=3/5` by default): the solvent favors
  chain-solvent contacts, swelling the chain relative to ideal. This is
  Flory's original 1953 mean-field estimate; a more precise
  renormalization-group value, :math:`\nu\approx0.588` (de Gennes,
  1979), is available separately via
  :meth:`RealChain.good_solvent_renormalization_group` for comparison --
  the two differ by about 2%.
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

import numpy as np

from chemistrykit.polymer.core.base_system import PolymerChainModel

__all__ = [
    "FLORY_EXPONENTS",
    "FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP",
    "flory_exponent",
    "IdealChain",
    "RealChain",
    "freely_jointed_chain",
    "worm_like_chain_mean_square_end_to_end",
]

#: dict: Standard Flory exponents nu for R ~ b*n**nu, by solvent quality
#: (Rubinstein & Colby, *Polymer Physics*, 2003, Table 3.2). "good" is
#: Flory's original 1953 mean-field estimate; see
#: FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP for the more precise
#: renormalization-group value.
FLORY_EXPONENTS = {
    "theta": 0.5,
    "good": 3.0 / 5.0,
    "poor": 1.0 / 3.0,
}

#: float: The renormalization-group-refined Flory exponent for a good
#: solvent in 3D, :math:`\nu\approx0.588` (P.-G. de Gennes, "Exponents
#: for the Excluded Volume Problem as Derived by the Wilson Method,"
#: Phys. Lett. A 38 (1972), 339-340; J. C. Le Guillou & J. Zinn-Justin,
#: Phys. Rev. Lett. 39 (1977), 95), refining Flory's original mean-field
#: estimate of exactly 3/5 (FLORY_EXPONENTS["good"]) by less than 2%.
FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP = 0.588


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
        """Build a :class:`RealChain` with Flory's mean-field good-solvent exponent (nu=3/5).

        See :meth:`good_solvent_renormalization_group` for the more
        precise renormalization-group value.
        """
        return cls(FLORY_EXPONENTS["good"])

    @classmethod
    def good_solvent_renormalization_group(cls) -> RealChain:
        r"""Build a :class:`RealChain` with the renormalization-group good-solvent exponent (nu ~= 0.588).

        This is de Gennes' (1979) more precise value, refining Flory's
        original mean-field estimate (:meth:`good_solvent`, nu=3/5) by
        about 2% -- see the module docstring.

        Examples
        --------
        >>> flory = RealChain.good_solvent()
        >>> de_gennes = RealChain.good_solvent_renormalization_group()
        >>> relative_difference = abs(flory.nu - de_gennes.nu) / flory.nu
        >>> round(float(relative_difference) * 100, 2)
        2.0
        """
        return cls(FLORY_EXPONENT_GOOD_SOLVENT_RENORMALIZATION_GROUP)

    @classmethod
    def poor_solvent(cls) -> RealChain:
        """Build a :class:`RealChain` with the poor-solvent Flory exponent (nu=1/3)."""
        return cls(FLORY_EXPONENTS["poor"])

    def mean_square_end_to_end(self, n, b):
        return (b * n**self.nu) ** 2

    def mean_square_radius_of_gyration(self, n, b):
        """See the module docstring: the 1/6 prefactor is an ideal-chain approximation, not exact for real chains."""
        return self.mean_square_end_to_end(n, b) / 6.0


def freely_jointed_chain(n: int, b: float, n_chains: int = 1, rng=None) -> np.ndarray:
    r"""Sample 3D conformations of Kuhn's freely jointed chain: `n` bonds of length `b` in uniformly random directions.

    Each bond vector is drawn independently and isotropically on the
    sphere of radius `b` (W. Kuhn, *Kolloid-Z.* 68, 2 (1934)); the chain
    starts at the origin. Averaged over many samples, the squared
    end-to-end distance converges to the exact ideal-chain result
    :math:`\langle R^2\rangle=nb^2` (:class:`IdealChain`), since the
    cross terms :math:`\langle\mathbf{b}_i\cdot\mathbf{b}_j\rangle`
    vanish for independent bonds.

    Parameters
    ----------
    n : int
        Number of bonds (Kuhn segments).
    b : float
        Bond (Kuhn) length.
    n_chains : int, optional
        Number of independent conformations to sample.
    rng : int, numpy.random.Generator, or None, optional
        Seed or generator for reproducibility.

    Returns
    -------
    ndarray, shape (n_chains, n + 1, 3)
        Bead positions of every sampled chain.

    Examples
    --------
    Every bond has exactly length `b`:

    >>> X = freely_jointed_chain(50, 0.5, n_chains=3, rng=0)
    >>> X.shape
    (3, 51, 3)
    >>> bool(np.allclose(np.linalg.norm(np.diff(X, axis=1), axis=-1), 0.5))
    True

    The sample mean of :math:`R^2` approaches :math:`nb^2=100\cdot1^2`:

    >>> X = freely_jointed_chain(100, 1.0, n_chains=20000, rng=1)
    >>> R2 = np.sum(X[:, -1] ** 2, axis=-1)
    >>> bool(abs(R2.mean() / 100.0 - 1.0) < 0.03)
    True
    """
    rng = np.random.default_rng(rng)
    bonds = rng.normal(size=(n_chains, n, 3))
    bonds *= b / np.linalg.norm(bonds, axis=-1, keepdims=True)
    positions = np.zeros((n_chains, n + 1, 3))
    positions[:, 1:] = np.cumsum(bonds, axis=1)
    return positions


def worm_like_chain_mean_square_end_to_end(L, P):
    r"""Kratky-Porod worm-like chain mean-square end-to-end distance.

    For a semiflexible chain of contour length :math:`L` whose tangent
    correlations decay as :math:`e^{-s/P}` along the contour (persistence
    length :math:`P`), Kratky and Porod (*Recl. Trav. Chim. Pays-Bas* 68,
    1106 (1949)) obtained

    .. math::

        \langle R^2\rangle = 2PL\left[1-\frac{P}{L}\left(1-e^{-L/P}\right)\right]

    which interpolates between a rigid rod, :math:`\langle R^2\rangle\to
    L^2` for :math:`L\ll P`, and an ideal random coil,
    :math:`\langle R^2\rangle\to2PL` for :math:`L\gg P` -- i.e. a
    freely jointed chain with Kuhn length :math:`b=2P`.

    Parameters
    ----------
    L : float or array-like of float
        Contour length (> 0).
    P : float
        Persistence length (> 0).

    Returns
    -------
    float or ndarray

    Examples
    --------
    A very short chain is a rigid rod, :math:`R^2\approx L^2`:

    >>> round(float(worm_like_chain_mean_square_end_to_end(1e-3, 50.0) / 1e-6), 4)
    1.0

    A very long chain is an ideal coil with Kuhn length :math:`2P`:

    >>> round(float(worm_like_chain_mean_square_end_to_end(1e6, 50.0) / (2 * 50.0 * 1e6)), 4)
    1.0
    """
    L = np.asarray(L, dtype=float)
    x = L / P
    return 2.0 * P * L * (1.0 - (1.0 - np.exp(-x)) / x)
