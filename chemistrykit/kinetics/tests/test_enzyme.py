"""Tests for chemistrykit.kinetics.systems.enzyme."""

import numpy as np
import pytest

from chemistrykit.kinetics.systems.enzyme import (
    MichaelisMentenProgress,
    competitive_inhibition_rate,
    fit_lineweaver_burk,
    michaelis_menten_rate,
    noncompetitive_inhibition_rate,
)


def test_michaelis_menten_rate_at_Km_is_half_Vmax():
    """By definition, Km is the substrate concentration giving v = Vmax/2."""
    assert michaelis_menten_rate(S=3.0, Vmax=12.0, Km=3.0) == pytest.approx(6.0)


def test_michaelis_menten_rate_saturates_to_Vmax():
    v = michaelis_menten_rate(S=1e6, Vmax=10.0, Km=1.0)
    assert v == pytest.approx(10.0, rel=1e-4)


def test_michaelis_menten_rate_is_first_order_at_low_substrate():
    """v ~ (Vmax/Km)*S when S << Km (first-order limit)."""
    Vmax, Km = 10.0, 100.0
    S = 0.01
    v = michaelis_menten_rate(S=S, Vmax=Vmax, Km=Km)
    assert v == pytest.approx((Vmax / Km) * S, rel=1e-3)


def test_fit_lineweaver_burk_recovers_exact_parameters():
    Vmax_true, Km_true = 15.0, 4.0
    S = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
    v = michaelis_menten_rate(S, Vmax_true, Km_true)
    fit = fit_lineweaver_burk(S, v)
    assert fit.Vmax == pytest.approx(Vmax_true, rel=1e-6)
    assert fit.Km == pytest.approx(Km_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)


def test_competitive_inhibition_raises_apparent_Km_but_not_Vmax():
    Vmax, Km, Ki = 10.0, 2.0, 1.0
    S = np.array([0.5, 1.0, 5.0, 50.0, 500.0])
    v_uninhibited = michaelis_menten_rate(S, Vmax, Km)
    v_inhibited = competitive_inhibition_rate(S, I=2.0, Vmax=Vmax, Km=Km, Ki=Ki)
    # At low [S] competitive inhibition suppresses the rate...
    assert v_inhibited[0] < v_uninhibited[0]
    # ...but at saturating [S] both approach the *same* Vmax (competition
    # can be out-competed by enough substrate).
    assert v_inhibited[-1] == pytest.approx(v_uninhibited[-1], rel=1e-2)


def test_noncompetitive_inhibition_lowers_apparent_Vmax():
    Vmax, Km, Ki = 10.0, 2.0, 1.0
    S = np.array([1e5])
    v_uninhibited = michaelis_menten_rate(S, Vmax, Km)
    v_inhibited = noncompetitive_inhibition_rate(S, I=1.0, Vmax=Vmax, Km=Km, Ki=Ki)
    # I == Ki halves the apparent Vmax, and at saturating S, v -> Vmax_app.
    assert v_inhibited[0] == pytest.approx(v_uninhibited[0] / 2.0, rel=1e-3)


def test_michaelis_menten_progress_reduces_to_first_order_at_low_substrate():
    """S0 << Km: d[S]/dt ~ -(Vmax/Km)*S, i.e. first-order decay."""
    Vmax, Km, S0 = 1.0, 1000.0, 1.0
    k_eff = Vmax / Km
    network = MichaelisMentenProgress(S0=S0, Vmax=Vmax, Km=Km)
    result = network.integrate((0.0, 5.0), dt=1e-3, method="rk4")
    expected = S0 * np.exp(-k_eff * result.t)
    np.testing.assert_allclose(result.concentration("S"), expected, atol=1e-4)


def test_michaelis_menten_progress_reduces_to_zero_order_at_high_substrate():
    """S0 >> Km: d[S]/dt ~ -Vmax, i.e. zero-order decay, until S approaches Km."""
    Vmax, Km, S0 = 1.0, 0.01, 10.0
    network = MichaelisMentenProgress(S0=S0, Vmax=Vmax, Km=Km)
    result = network.integrate((0.0, 1.0), dt=1e-4, method="rk4")
    # Over a short early window S is still >> Km, so depletion is
    # essentially linear at rate Vmax.
    early = result.t < 1.0
    expected_early = S0 - Vmax * result.t[early]
    np.testing.assert_allclose(result.concentration("S")[early], expected_early, atol=5e-3)


def test_michaelis_menten_progress_never_goes_negative():
    network = MichaelisMentenProgress(S0=1.0, Vmax=5.0, Km=0.1)
    result = network.integrate((0.0, 5.0), dt=1e-4, method="rk4")
    assert np.all(result.concentration("S") >= -1e-6)
