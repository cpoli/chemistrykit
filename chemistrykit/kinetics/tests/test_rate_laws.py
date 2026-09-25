"""Tests for chemistrykit.kinetics.systems.rate_laws against closed-form results."""

import numpy as np
import pytest

from chemistrykit.kinetics.systems.rate_laws import FirstOrder, SecondOrder, ZeroOrder


def test_zero_order_concentration_and_half_life():
    law = ZeroOrder(k=0.1, C0=1.0)
    assert law.concentration(0.0) == pytest.approx(1.0)
    assert law.concentration(5.0) == pytest.approx(0.5)
    assert law.half_life() == pytest.approx(5.0)
    # By definition, concentration at the half-life is exactly half C0.
    assert law.concentration(law.half_life()) == pytest.approx(law.C0 / 2.0)


def test_zero_order_rate_is_constant_until_exhaustion():
    law = ZeroOrder(k=0.3, C0=2.0)  # exhausted at t = C0/k = 6.67
    assert law.rate() == pytest.approx(0.3)
    assert law.rate(t=5.0) == pytest.approx(0.3)
    assert law.rate(t=100.0) == 0.0
    np.testing.assert_allclose(law.rate(np.array([0.0, 6.0, 7.0])), [0.3, 0.3, 0.0])


def test_zero_order_concentration_never_goes_negative():
    law = ZeroOrder(k=0.1, C0=1.0)
    assert law.concentration(10.0) == pytest.approx(0.0)
    assert law.concentration(50.0) == 0.0
    np.testing.assert_allclose(law.concentration(np.array([0.0, 5.0, 15.0])), [1.0, 0.5, 0.0])


def test_first_order_matches_exponential_decay():
    k, C0 = 0.7, 3.0
    law = FirstOrder(k=k, C0=C0)
    t = np.linspace(0.0, 10.0, 50)
    np.testing.assert_allclose(law.concentration(t), C0 * np.exp(-k * t))


def test_first_order_half_life_independent_of_initial_concentration():
    """The defining signature of first-order kinetics: t_1/2 = ln(2)/k regardless of C0."""
    k = 0.5
    t_half_1 = FirstOrder(k=k, C0=1.0).half_life()
    t_half_2 = FirstOrder(k=k, C0=100.0).half_life()
    assert t_half_1 == pytest.approx(t_half_2)
    assert t_half_1 == pytest.approx(np.log(2.0) / k)


def test_first_order_concentration_at_half_life_is_half_c0():
    law = FirstOrder(k=1.3, C0=5.0)
    assert law.concentration(law.half_life()) == pytest.approx(law.C0 / 2.0, rel=1e-9)


def test_second_order_matches_integrated_law():
    k, C0 = 0.4, 2.0
    law = SecondOrder(k=k, C0=C0)
    t = np.linspace(0.0, 10.0, 50)
    expected = C0 / (1.0 + k * C0 * t)
    np.testing.assert_allclose(law.concentration(t), expected)
    # Equivalent form: 1/[A] = 1/[A]_0 + kt.
    np.testing.assert_allclose(1.0 / law.concentration(t), 1.0 / C0 + k * t)


def test_second_order_half_life_depends_on_initial_concentration():
    """Unlike first order, t_1/2 = 1/(k*C0) shrinks as C0 grows."""
    k = 0.4
    t_half_low = SecondOrder(k=k, C0=1.0).half_life()
    t_half_high = SecondOrder(k=k, C0=10.0).half_life()
    assert t_half_high < t_half_low
    assert t_half_low == pytest.approx(1.0 / (k * 1.0))


@pytest.mark.parametrize("cls", [ZeroOrder, FirstOrder, SecondOrder])
def test_invalid_parameters_are_rejected(cls):
    with pytest.raises(ValueError):
        cls(k=-1.0, C0=1.0)
    with pytest.raises(ValueError):
        cls(k=1.0, C0=-1.0)
