"""Tests for chemistrykit.statmech.systems.maxwell_boltzmann against closed-form results."""

import numpy as np
import pytest
from scipy.integrate import quad

from chemistrykit.statmech.systems.maxwell_boltzmann import (
    MaxwellBoltzmannSpeedDistribution,
    mean_speed,
    most_probable_speed,
    rms_speed,
)

MASS = 6.63e-26  # ~argon, kg
TEMPERATURE = 298.15


def test_characteristic_speed_ratios_are_exact():
    """v_p : v_bar : v_rms = sqrt(2) : sqrt(8/pi) : sqrt(3), independent of mass/T."""
    vp = most_probable_speed(MASS, TEMPERATURE)
    vbar = mean_speed(MASS, TEMPERATURE)
    vrms = rms_speed(MASS, TEMPERATURE)
    np.testing.assert_allclose(vbar / vp, np.sqrt(8.0 / np.pi) / np.sqrt(2.0), rtol=1e-10)
    np.testing.assert_allclose(vrms / vp, np.sqrt(3.0) / np.sqrt(2.0), rtol=1e-10)


def test_speed_ratios_independent_of_mass_and_temperature():
    ratio_1 = mean_speed(MASS, TEMPERATURE) / most_probable_speed(MASS, TEMPERATURE)
    ratio_2 = mean_speed(2.0 * MASS, 500.0) / most_probable_speed(2.0 * MASS, 500.0)
    assert ratio_1 == pytest.approx(ratio_2)


def test_pdf_integrates_to_one():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    total, _ = quad(dist.pdf, 0.0, np.inf)
    assert total == pytest.approx(1.0, abs=1e-6)


def test_pdf_peaks_at_most_probable_speed():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    vp = dist.most_probable_speed()
    v = np.linspace(0.1, 5.0 * vp, 10000)
    v_at_peak = v[np.argmax(dist.pdf(v))]
    assert v_at_peak == pytest.approx(vp, rel=2e-3)


def test_cdf_zero_at_origin_and_one_at_infinity():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    assert float(dist.cdf(0.0)) == pytest.approx(0.0, abs=1e-12)
    assert float(dist.cdf(1.0e5)) == pytest.approx(1.0, abs=1e-6)


def test_cdf_derivative_matches_pdf():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    v = 400.0
    h = 1e-2
    numerical_derivative = (float(dist.cdf(v + h)) - float(dist.cdf(v - h))) / (2.0 * h)
    assert numerical_derivative == pytest.approx(float(dist.pdf(v)), rel=1e-3)


def test_sample_mean_matches_analytic_mean_speed():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    speeds = dist.sample(200_000, rng=0)
    assert speeds.mean() == pytest.approx(dist.mean_speed(), rel=0.01)


def test_sample_rms_matches_analytic_rms_speed():
    dist = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=TEMPERATURE)
    speeds = dist.sample(200_000, rng=1)
    rms = np.sqrt(np.mean(speeds**2))
    assert rms == pytest.approx(dist.rms_speed(), rel=0.01)


def test_higher_temperature_gives_higher_characteristic_speeds():
    dist_cold = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=100.0)
    dist_hot = MaxwellBoltzmannSpeedDistribution(mass=MASS, temperature=1000.0)
    assert dist_hot.mean_speed() > dist_cold.mean_speed()


def test_rejects_nonpositive_mass_or_temperature():
    with pytest.raises(ValueError):
        MaxwellBoltzmannSpeedDistribution(mass=0.0, temperature=300.0)
    with pytest.raises(ValueError):
        MaxwellBoltzmannSpeedDistribution(mass=1.0, temperature=-1.0)


def test_md_lj_fluid_speeds_match_maxwell_boltzmann_at_equilibrium():
    """Cross-domain check: an equilibrated LJFluid's speed histogram follows this distribution."""
    from chemistrykit.md.systems.lj_fluid import LJFluid

    fluid = LJFluid.from_lattice(n_per_side=6, density=0.4, temperature=1.5, rng=42)
    fluid.run(dt=0.002, n_steps=500, sample_every=500)  # brief equilibration
    speeds = np.linalg.norm(fluid.velocities, axis=1)
    dist = MaxwellBoltzmannSpeedDistribution(mass=1.0, temperature=fluid.temperature(), k_b=1.0)
    assert speeds.mean() == pytest.approx(dist.mean_speed(), rel=0.15)
