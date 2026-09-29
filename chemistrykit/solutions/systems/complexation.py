r"""Stepwise metal-ligand complexation: species distribution, Bjerrum's formation function, and mass balance.

A metal ion M binding up to `N` ligands L forms :math:`ML_n` with stepwise
constants :math:`K_n = [ML_n]/([ML_{n-1}][L])` and cumulative constants
:math:`\beta_n = K_1K_2\cdots K_n = [ML_n]/([M][L]^n)`. The fraction of
total metal present as :math:`ML_n` depends only on the free-ligand
concentration,

.. math::

    \alpha_n = \frac{\beta_n[L]^n}{\sum_{k=0}^N\beta_k[L]^k}, \qquad \beta_0 = 1

-- the complexation analogue of the polyprotic-acid fractions in
:func:`chemistrykit.solutions.polyprotic_fractions`, with :math:`[L]` in
place of :math:`1/[H^+]`. The mean number of bound ligands,
:math:`\bar n = \sum_n n\alpha_n`, is Bjerrum's formation function (J.
Bjerrum, *Metal Ammine Formation in Aqueous Solution*, Copenhagen: Haase,
1941; Harris, *Quantitative Chemical Analysis*, 9th ed., Ch. 6-4 and
12-1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "cumulative_formation_constants",
    "complex_fractions",
    "average_ligand_number",
    "ComplexationEquilibrium",
    "solve_complexation",
]


def cumulative_formation_constants(K_stepwise) -> np.ndarray:
    r"""Cumulative constants :math:`\beta_n=\prod_{i\le n}K_i` from stepwise ones.

    Parameters
    ----------
    K_stepwise : array-like of float
        :math:`K_1, \dots, K_N`.

    Returns
    -------
    ndarray, shape (N,)

    Examples
    --------
    >>> cumulative_formation_constants([100.0, 10.0]).tolist()
    [100.0, 1000.0]
    """
    return np.cumprod(np.asarray(K_stepwise, dtype=np.float64))


def complex_fractions(L_free, beta) -> np.ndarray:
    r"""Fractions :math:`\alpha_0, \dots, \alpha_N` of total metal present as :math:`M, ML, \dots, ML_N`.

    Parameters
    ----------
    L_free : float or array-like of float
        Free (unbound) ligand concentration(s), in mol/L.
    beta : array-like of float
        Cumulative formation constants :math:`\beta_1, \dots, \beta_N`.

    Returns
    -------
    ndarray
        Shape ``(N + 1, len(L_free))``; row ``n`` is :math:`\alpha_n`. The
        rows sum to 1.

    Examples
    --------
    Ag+ with NH3 (:math:`\log\beta_1=3.31`, :math:`\log\beta_2=7.23`):
    at :math:`[NH_3]=0.01` M almost all silver is the diammine complex:

    >>> alpha = complex_fractions(0.01, [10**3.31, 10**7.23])
    >>> [round(float(a), 3) for a in alpha[:, 0]]
    [0.001, 0.012, 0.988]
    """
    L = np.atleast_1d(np.asarray(L_free, dtype=np.float64))
    beta_all = np.concatenate(([1.0], np.asarray(beta, dtype=np.float64)))
    terms = np.array([b * L**n for n, b in enumerate(beta_all)])
    return terms / terms.sum(axis=0)


def average_ligand_number(L_free, beta):
    r"""Bjerrum's formation function :math:`\bar n=\sum_n n\alpha_n`, the mean number of ligands per metal.

    Parameters
    ----------
    L_free : float or array-like of float
        Free ligand concentration(s), in mol/L.
    beta : array-like of float
        Cumulative formation constants.

    Returns
    -------
    ndarray, shape (len(L_free),)

    Examples
    --------
    For a single step, :math:`\bar n = 1/2` exactly where :math:`[L]=1/K_1`:

    >>> round(float(average_ligand_number(1e-3, [1e3])[0]), 6)
    0.5
    """
    alpha = complex_fractions(L_free, beta)
    n = np.arange(alpha.shape[0])[:, None]
    return (n * alpha).sum(axis=0)


@dataclass
class ComplexationEquilibrium:
    """Result of :func:`solve_complexation`."""

    free_ligand: float
    """float: Free ligand concentration [L], in mol/L."""

    species: np.ndarray
    """ndarray, shape (N + 1,): Concentrations of M, ML, ..., ML_N, in mol/L."""

    average_ligand_number: float
    """float: Bjerrum's formation function at equilibrium."""

    @property
    def free_metal(self) -> float:
        """float: Free metal concentration [M], in mol/L."""
        return float(self.species[0])


def solve_complexation(M_total: float, L_total: float, beta) -> ComplexationEquilibrium:
    r"""Equilibrium speciation from total metal and total ligand concentrations.

    Solves the ligand mass balance :math:`L_T = [L] + M_T\,\bar n([L])` for
    the free ligand concentration by bracketed root-finding
    (:func:`scipy.optimize.brentq` on :math:`0 < [L] \le L_T`). The left
    side minus the right is monotonic in :math:`[L]`, so the root is
    unique. Protonation of the ligand and hydrolysis of the metal are
    neglected; fold them in through conditional constants if needed.

    Parameters
    ----------
    M_total : float
        Total metal concentration, in mol/L.
    L_total : float
        Total ligand concentration, in mol/L.
    beta : array-like of float
        Cumulative formation constants.

    Returns
    -------
    ComplexationEquilibrium

    Examples
    --------
    0.001 M Ag+ in 0.1 M NH3: nearly all silver ends up as Ag(NH3)2+:

    >>> eq = solve_complexation(1e-3, 0.1, [10**3.31, 10**7.23])
    >>> round(eq.average_ligand_number, 3)
    1.999
    """
    if M_total < 0 or L_total <= 0:
        raise ValueError("M_total must be non-negative and L_total positive")
    beta = np.asarray(beta, dtype=np.float64)

    def residual(L):
        return L + M_total * float(average_ligand_number(L, beta)[0]) - L_total

    L_free = brentq(residual, 0.0, L_total, xtol=1e-15 * max(L_total, 1.0), rtol=1e-12)
    alpha = complex_fractions(L_free, beta)[:, 0]
    return ComplexationEquilibrium(free_ligand=float(L_free), species=M_total * alpha, average_ligand_number=float(np.dot(np.arange(alpha.size), alpha)))
