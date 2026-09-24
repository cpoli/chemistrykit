"""Tests for chemistrykit.surface.systems.eley_rideal against closed-form results."""

import numpy as np
import pytest

from chemistrykit.surface.systems.eley_rideal import er_rate
from chemistrykit.surface.systems.langmuir_hinshelwood import lh_rate_dual_site


def test_er_rate_closed_form():
    k, K_A, P_A, P_B = 3.0, 2.0, 0.7, 1.5
    assert er_rate(k, K_A, P_A, P_B) == pytest.approx(k * K_A * P_A * P_B / (1 + K_A * P_A))


def test_er_rate_is_monotonic_in_P_A_and_saturates_at_k_P_B():
    P_A = np.linspace(0.01, 50.0, 200)
    rate = er_rate(4.0, 2.0, P_A, 0.2)
    assert np.all(np.diff(rate) > 0)
    assert er_rate(4.0, 2.0, 1e9, 0.2) == pytest.approx(0.8, rel=1e-6)


def test_er_rate_is_first_order_in_P_B_unlike_lh():
    P_B = np.array([0.1, 1.0, 10.0])
    np.testing.assert_allclose(er_rate(4.0, 2.0, 1.0, P_B) / P_B, er_rate(4.0, 2.0, 1.0, 1.0))
    lh = lh_rate_dual_site(4.0, 2.0, 1.0, 5.0, P_B)
    assert not np.allclose(lh / P_B, lh[0] / P_B[0])
