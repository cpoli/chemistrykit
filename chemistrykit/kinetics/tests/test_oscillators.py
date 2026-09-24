"""Tests for chemistrykit.kinetics.systems.oscillators.Brusselator."""

import numpy as np

from chemistrykit.kinetics.systems.oscillators import Brusselator


def test_fixed_point_is_a_zero_of_the_vector_field():
    for A, B in [(1.0, 3.0), (2.0, 5.0), (0.5, 1.0)]:
        system = Brusselator(A=A, B=B)
        fp = system.fixed_point()
        np.testing.assert_allclose(system.rhs(fp, 0.0), 0.0, atol=1e-12)


def test_hopf_threshold_matches_definition():
    assert Brusselator(A=1.0, B=3.0).is_above_hopf_threshold() is True
    assert Brusselator(A=1.0, B=1.5).is_above_hopf_threshold() is False
    assert Brusselator(A=2.0, B=5.0).is_above_hopf_threshold() is False  # exactly at threshold: 5 == 1 + 2^2
    assert Brusselator(A=2.0, B=5.01).is_above_hopf_threshold() is True


def test_below_hopf_threshold_relaxes_to_fixed_point():
    """Below the Hopf threshold the fixed point is a stable focus: starting
    away from it, the trajectory should settle back down."""
    system = Brusselator(X0=2.0, Y0=0.5, A=1.0, B=1.2)
    fp = system.fixed_point()
    result = system.integrate((0.0, 200.0), dt=1e-2, method="rk4")
    final = result.y[-1]
    np.testing.assert_allclose(final, fp, atol=5e-2)


def test_above_hopf_threshold_sustains_oscillation():
    """Above the Hopf threshold, the fixed point is unstable and the
    trajectory should settle into a limit cycle: concentrations keep
    oscillating with an amplitude that does not decay to zero."""
    system = Brusselator(X0=1.0, Y0=1.0, A=1.0, B=3.0)
    result = system.integrate((0.0, 200.0), dt=1e-2, method="rk4")
    # Look only at the second half of the run, well past any transient.
    tail = result.y[len(result.y) // 2 :]
    X_tail = tail[:, 0]
    amplitude = np.max(X_tail) - np.min(X_tail)
    assert amplitude > 0.5  # sustained oscillation, not decayed to the fixed point


def test_concentrations_stay_non_negative_above_threshold():
    system = Brusselator(X0=0.1, Y0=0.1, A=1.0, B=3.0)
    result = system.integrate((0.0, 100.0), dt=1e-2, method="rk4")
    assert np.all(result.y >= -1e-6)


def test_oregonator_fixed_point_is_zero_of_vector_field():
    from chemistrykit.kinetics.systems.oscillators import Oregonator

    for f in (0.5, 1.0, 3.0):
        system = Oregonator(f=f)
        np.testing.assert_allclose(system.rhs(system.fixed_point(), 0.0), 0.0, atol=1e-8)


def test_oregonator_oscillates_for_f1_and_is_stable_for_f3():
    from chemistrykit.kinetics.systems.oscillators import Oregonator

    osc = Oregonator(f=1.0).integrate((0.0, 60.0), method="dopri5", max_steps=2_000_000)
    x_late = osc.y[osc.t > 30.0, 0]
    assert x_late.max() - x_late.min() > 0.5

    stable = Oregonator(f=3.0)
    result = stable.integrate((0.0, 60.0), method="dopri5", max_steps=2_000_000)
    np.testing.assert_allclose(result.y[-1], stable.fixed_point(), rtol=1e-3)
