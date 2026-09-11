r"""Exact combinatorial helpers for canonical-ensemble microstate counting.

Used by :class:`chemistrykit.statmech.systems.lattice_gas.LatticeGasAdsorption`
to count the number of ways to arrange `N` indistinguishable adsorbed
molecules among `M` distinguishable lattice sites.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

__all__ = ["ln_factorial", "ln_binomial"]


def ln_factorial(n):
    r"""Return :math:`\ln(n!)`, via the exact log-gamma function :math:`\ln\Gamma(n+1)`.

    Traditional statistical-mechanics derivations of combinatorial
    entropy approximate :math:`\ln n!` with Stirling's approximation
    (:math:`\ln n! \approx n\ln n - n`); chemistrykit instead uses
    SciPy's numerically exact ``gammaln`` (accurate to floating-point
    precision for every :math:`n\geq0`, including small `n` where
    Stirling's approximation is poor), so no accuracy is sacrificed for
    the sake of matching the textbook derivation -- the two agree closely
    only in the large-`n` limit that :func:`ln_binomial`'s docstring
    checks explicitly.

    Parameters
    ----------
    n : float or array-like of float
        Need not be an integer (:math:`\Gamma(n+1)` generalizes `n!`).

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> import numpy as np
    >>> round(float(ln_factorial(5)), 6) == round(float(np.log(120)), 6)
    True
    """
    return gammaln(np.asarray(n, dtype=np.float64) + 1.0)


def ln_binomial(n, k):
    r"""Return :math:`\ln\binom{n}{k} = \ln n! - \ln k! - \ln(n-k)!`.

    Parameters
    ----------
    n, k : float or array-like of float

    Returns
    -------
    float or ndarray

    Examples
    --------
    Agrees with direct combinatorics for small, exact integers:

    >>> import numpy as np
    >>> round(float(np.exp(ln_binomial(5, 2))), 6)
    10.0

    And, for large `n`, the peak degeneracy at :math:`k=n/2` approaches
    Stirling's classic large-`n` estimate :math:`\ln\binom{n}{n/2}\approx
    n\ln2` (the basis of the large-`M` limit checked in
    :mod:`chemistrykit.statmech.systems.lattice_gas`):

    >>> n = 100_000
    >>> relative_error = abs(ln_binomial(n, n // 2) - n * np.log(2)) / (n * np.log(2))
    >>> bool(relative_error < 1e-4)
    True
    """
    n = np.asarray(n, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    return ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
