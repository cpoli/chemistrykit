"""Tests for chemistrykit.kinetics.systems.rate_theory against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import K_B, NA, H, R
from chemistrykit.kinetics.systems.arrhenius import fit_arrhenius
from chemistrykit.kinetics.systems.rate_theory import (
    collision_theory_rate_constant,
    diffusion_limited_rate_constant,
    eyring_rate_constant,
    fit_eyring,
    smoluchowski_rate_constant,
    smoluchowski_transient_rate_constant,
)


def test_collision_theory_matches_formula():
    T, sigma, mu, Ea, P = 500.0, 3e-19, 2e-26, 40e3, 0.2
    expected = P * sigma * np.sqrt(8 * K_B * T / (np.pi * mu)) * NA * np.exp(-Ea / (R * T))
    assert collision_theory_rate_constant(T, sigma, mu, Ea=Ea, steric_factor=P) == pytest.approx(expected, rel=1e-12)


def test_collision_theory_apparent_activation_energy_is_Ea_plus_half_RT():
    """d ln k / d(1/T) = -(Ea + RT/2)/R for k ∝ sqrt(T) exp(-Ea/RT)."""
    Ea, T = 50e3, 400.0
    h = 1e-3
    k_hi = collision_theory_rate_constant(T + h, 4e-19, 1e-26, Ea=Ea)
    k_lo = collision_theory_rate_constant(T - h, 4e-19, 1e-26, Ea=Ea)
    slope = (np.log(k_hi) - np.log(k_lo)) / (1 / (T + h) - 1 / (T - h))
    assert -slope * R == pytest.approx(Ea + 0.5 * R * T, rel=1e-6)


def test_smoluchowski_formula_and_transient_limit():
    D, Rc = 3e-9, 4e-10
    kD = smoluchowski_rate_constant(D, Rc)
    assert kD == pytest.approx(4 * np.pi * D * Rc * NA)
    assert smoluchowski_transient_rate_constant(D, Rc, 1e3) == pytest.approx(kD, rel=1e-5)
    assert smoluchowski_transient_rate_constant(D, Rc, 1e-12) > kD


def test_diffusion_limited_equals_smoluchowski_with_stokes_einstein():
    T, eta, r = 298.15, 8.9e-4, 3e-10
    D_each = K_B * T / (6 * np.pi * eta * r)
    assert diffusion_limited_rate_constant(T, eta) == pytest.approx(smoluchowski_rate_constant(2 * D_each, 2 * r), rel=1e-12)


def test_eyring_zero_barrier_is_universal_frequency():
    assert eyring_rate_constant(310.0, 0.0, 0.0) == pytest.approx(K_B * 310.0 / H)


def test_fit_eyring_recovers_parameters_from_noisy_data():
    rng = np.random.default_rng(3)
    T = np.linspace(290.0, 350.0, 12)
    k = eyring_rate_constant(T, dH=75e3, dS=-30.0) * (1 + rng.normal(0, 0.01, T.shape))
    fit = fit_eyring(T, k)
    assert fit.dH == pytest.approx(75e3, rel=0.02)
    assert fit.dS == pytest.approx(-30.0, abs=5.0)
    assert fit.dG(300.0) == pytest.approx(fit.dH - 300.0 * fit.dS)


def test_arrhenius_activation_energy_is_dH_plus_RT():
    """For Eyring kinetics, Ea = dH + R*T_mid (to first order in the T range)."""
    T = np.linspace(295.0, 305.0, 11)
    dH = 60e3
    fit = fit_arrhenius(T, eyring_rate_constant(T, dH=dH, dS=10.0))
    assert fit.Ea == pytest.approx(dH + R * 300.0, rel=1e-4)
