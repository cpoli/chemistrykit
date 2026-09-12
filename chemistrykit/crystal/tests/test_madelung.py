"""Tests for chemistrykit.crystal.systems.madelung: convergence to the known NaCl value."""

import numpy as np
import pytest

from chemistrykit.crystal.systems.madelung import MADELUNG_CONSTANT_NACL_LITERATURE, madelung_constant_nacl
from chemistrykit.crystal.utils.lattice_sums import evjen_lattice_sum_cubic_alternating


def test_madelung_constant_nacl_matches_literature_value():
    assert madelung_constant_nacl(n_shells=15) == pytest.approx(MADELUNG_CONSTANT_NACL_LITERATURE, abs=1e-4)


def test_madelung_constant_converges_as_n_shells_grows():
    values = [madelung_constant_nacl(n) for n in (2, 5, 10, 15, 20)]
    errors = [abs(v - MADELUNG_CONSTANT_NACL_LITERATURE) for v in values]
    # Not monotonic term-by-term, but the envelope of the error shrinks:
    # the largest error in the second half must be far smaller than the
    # largest error in the first half.
    assert max(errors[3:]) < max(errors[:2]) / 10.0


def test_evjen_weighted_sum_converges_much_faster_than_naive_truncation():
    def naive_sum(n):
        total = 0.0
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                for k in range(-n, n + 1):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    total += (-1) ** (i + j + k) / np.sqrt(i * i + j * j + k * k)
        return -total

    naive_values = [naive_sum(n) for n in range(5, 10)]
    evjen_values = [-evjen_lattice_sum_cubic_alternating(n) for n in range(5, 10)]
    assert np.std(evjen_values) < np.std(naive_values) / 100.0


def test_evjen_sum_rejects_nonpositive_n_shells():
    with pytest.raises(ValueError):
        evjen_lattice_sum_cubic_alternating(0)
