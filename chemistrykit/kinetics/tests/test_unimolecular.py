"""Tests for the Lindemann-Hinshelwood falloff model."""

import numpy as np
import pytest

from chemistrykit.kinetics.systems.unimolecular import LindemannHinshelwood


def test_limits_and_falloff_center():
    lh = LindemannHinshelwood(k1=2.0, k_minus1=5.0, k2=40.0)
    assert lh.rate_constant(1e-6) / 1e-6 == pytest.approx(lh.k0, rel=1e-5)
    assert lh.rate_constant(1e9) == pytest.approx(lh.k_inf, rel=1e-6)
    assert lh.rate_constant(lh.M_half) == pytest.approx(lh.k_inf / 2)


def test_reduced_pressure_form_and_lindemann_linearization():
    lh = LindemannHinshelwood(k1=2.0, k_minus1=5.0, k2=40.0)
    M = np.logspace(-2, 3, 11)
    Pr = lh.reduced_pressure(M)
    np.testing.assert_allclose(lh.rate_constant(M), lh.k_inf * Pr / (1 + Pr))
    np.testing.assert_allclose(lh.inverse_rate_constant(M), 1.0 / lh.rate_constant(M))


def test_full_mechanism_decays_at_steady_state_rate():
    lh = LindemannHinshelwood(k1=0.1, k_minus1=10.0, k2=100.0)
    M = 20.0
    result = lh.network(M).integrate((0.0, 20.0), dt=1e-3)
    A, t = result.concentration("A"), result.t
    i = len(t) // 2
    k_numeric = -np.log(A[-1] / A[i]) / (t[-1] - t[i])
    assert k_numeric == pytest.approx(lh.rate_constant(M), rel=5e-3)
    total = A + result.concentration("A*") + result.concentration("P")
    np.testing.assert_allclose(total, 1.0, atol=1e-9)
