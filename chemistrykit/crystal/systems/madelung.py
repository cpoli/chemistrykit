r"""The Madelung constant of the NaCl (rock-salt) structure, from a genuinely converging lattice sum.

The Madelung constant `M` of an ionic crystal is defined so that the
electrostatic lattice energy per ion pair is

.. math::

    U = -\frac{M z_+ z_- e^2}{4\pi\varepsilon_0 r_0}

with :math:`r_0` the nearest-neighbor (cation-anion) distance and
:math:`z_\pm` the (unsigned) ionic charge numbers (West, *Solid State
Chemistry and its Applications*, 2nd ed. (2014), Ch. 1.4; this is the `M`
consumed by :class:`chemistrykit.crystal.systems.lattice_energy.BornLandeLatticeEnergy`).
For NaCl, `M` is defined by the lattice sum

.. math::

    M = -\sum_{j\neq0}(-1)^{i+j+k}\,\frac{r_0}{r_j}

which, remarkably, can be evaluated over a *simple cubic* lattice of
alternating unit charges at spacing :math:`r_0`: NaCl's two
interpenetrating FCC sublattices (Na+ at the FCC points, Cl- at the FCC
points offset by half the conventional cell edge) combine into a single
simple cubic lattice of spacing :math:`a/2` (`a` the conventional cubic
cell edge) with ion charge alternating by the parity of :math:`i+j+k`,
exactly the object :func:`chemistrykit.crystal.utils.lattice_sums.evjen_lattice_sum_cubic_alternating`
sums (Ashcroft & Mermin, *Solid State Physics*, 1976, problem 20.1).

**Why this needs Evjen's method, not a naive truncated sum.** The series
above is only *conditionally* convergent -- it has no well-defined value
independent of summation order, so truncating the sum to a growing cube
or sphere of lattice points does not settle down to the right answer as
the cutoff grows; it *oscillates* (each newly-added shell is not
itself charge-neutral). Demonstrating this explicitly:

>>> import numpy as np
>>> def naive_sum(n):
...     total = 0.0
...     for i in range(-n, n + 1):
...         for j in range(-n, n + 1):
...             for k in range(-n, n + 1):
...                 if i == 0 and j == 0 and k == 0:
...                     continue
...                 total += (-1) ** (i + j + k) / np.sqrt(i * i + j * j + k * k)
...     return -total
>>> naive_values = [naive_sum(n) for n in range(6, 11)]
>>> bool(np.std(naive_values) > 0.05)  # still swinging wildly at n=6..10
True
>>> evjen_values = [madelung_constant_nacl(n) for n in range(6, 11)]
>>> bool(np.std(evjen_values) < 1e-4)  # Evjen's method has already settled down
True

whereas :func:`madelung_constant_nacl` (Evjen-weighted) converges quickly
and matches the literature value :math:`M_{NaCl}=1.747565` (Kittel,
*Introduction to Solid State Physics*, 8th ed., Table 1, Ch. 3) to 5
significant figures well before `n_shells=20`.
"""

from __future__ import annotations

from chemistrykit.crystal.utils.lattice_sums import evjen_lattice_sum_cubic_alternating

__all__ = ["MADELUNG_CONSTANT_NACL_LITERATURE", "madelung_constant_nacl"]

#: float: The accepted literature value of the NaCl Madelung constant
#: (Kittel, *Introduction to Solid State Physics*, 8th ed., Table 1, Ch.
#: 3), given here purely as a citation/regression-test reference -- see
#: the module docstring for the independent numerical derivation.
MADELUNG_CONSTANT_NACL_LITERATURE = 1.747565


def madelung_constant_nacl(n_shells: int = 12) -> float:
    r"""Compute the NaCl Madelung constant via Evjen's converging lattice summation.

    Parameters
    ----------
    n_shells : int, default 12
        Half-width (in simple-cubic lattice units) of the summation
        region passed to
        :func:`chemistrykit.crystal.utils.lattice_sums.evjen_lattice_sum_cubic_alternating`.
        The default already agrees with the literature value
        (:data:`MADELUNG_CONSTANT_NACL_LITERATURE`) to 5 significant
        figures.

    Returns
    -------
    float

    Examples
    --------
    >>> round(madelung_constant_nacl(), 4)
    1.7476
    >>> round(madelung_constant_nacl(20), 6)
    1.747565
    """
    return -evjen_lattice_sum_cubic_alternating(n_shells)
