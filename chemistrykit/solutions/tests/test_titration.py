"""Tests for chemistrykit.solutions.systems.titration against closed-form results."""

import numpy as np
import pytest

from chemistrykit.solutions.systems.titration import (
    StrongAcidStrongBaseTitration,
    WeakAcidStrongBaseTitration,
    WeakBaseStrongAcidTitration,
)


def test_strong_acid_strong_base_ph_at_start_matches_pure_strong_acid():
    titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    pH0 = titration.pH_at(np.array([0.0]))[0]
    assert pH0 == pytest.approx(1.0, abs=1e-3)


def test_strong_acid_strong_base_ph_at_equivalence_is_neutral():
    titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    Vb_eq = titration.equivalence_volume()
    pH_eq = titration.pH_at(np.array([Vb_eq]))[0]
    assert pH_eq == pytest.approx(7.0, abs=1e-6)


def test_strong_acid_strong_base_numeric_equivalence_matches_analytic():
    titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    Vb = np.linspace(1e-6, 0.09, 20000)
    Vb_numeric = titration.find_equivalence_point(Vb)
    assert Vb_numeric == pytest.approx(titration.equivalence_volume(), abs=1e-3)


def test_strong_acid_strong_base_ph_monotonically_increases():
    titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    Vb = np.linspace(0.0, 0.09, 200)
    pH = titration.pH_at(Vb)
    assert np.all(np.diff(pH) > 0)


def test_weak_acid_strong_base_half_equivalence_ph_equals_pka():
    Ka = 1.8e-5
    titration = WeakAcidStrongBaseTitration(Ca=0.100, Va=0.050, Ka=Ka, Cb=0.100)
    Vb_half = titration.equivalence_volume() / 2.0
    pH_half = titration.pH_at(np.array([Vb_half]))[0]
    assert pH_half == pytest.approx(-np.log10(Ka), abs=1e-3)


def test_weak_acid_strong_base_equivalence_point_is_basic():
    """Past the equivalence point, the solution is the conjugate base's salt: pH > 7."""
    titration = WeakAcidStrongBaseTitration(Ca=0.100, Va=0.050, Ka=1.8e-5, Cb=0.100)
    Vb_eq = titration.equivalence_volume()
    pH_eq = titration.pH_at(np.array([Vb_eq]))[0]
    assert pH_eq > 7.0


def test_weak_acid_strong_base_initial_ph_matches_weak_acid_alone():
    from chemistrykit.solutions.systems.acid_base import WeakAcid

    Ca, Ka = 0.100, 1.8e-5
    titration = WeakAcidStrongBaseTitration(Ca=Ca, Va=0.050, Ka=Ka, Cb=0.100)
    pH0 = titration.pH_at(np.array([0.0]))[0]
    expected = WeakAcid(Ca=Ca, Ka=Ka).pH()
    assert pH0 == pytest.approx(expected, abs=1e-3)


def test_weak_base_strong_acid_half_equivalence_poh_equals_pkb():
    Kb = 1.8e-5
    titration = WeakBaseStrongAcidTitration(Cb=0.100, Vb0=0.050, Kb=Kb, Ca=0.100)
    Va_half = titration.equivalence_volume() / 2.0
    pH_half = titration.pH_at(np.array([Va_half]))[0]
    pOH_half = 14.0 - pH_half
    assert pOH_half == pytest.approx(-np.log10(Kb), abs=1e-3)


def test_weak_base_strong_acid_equivalence_point_is_acidic():
    titration = WeakBaseStrongAcidTitration(Cb=0.100, Vb0=0.050, Kb=1.8e-5, Ca=0.100)
    Va_eq = titration.equivalence_volume()
    pH_eq = titration.pH_at(np.array([Va_eq]))[0]
    assert pH_eq < 7.0


def test_titration_curve_dataclass_shapes_match():
    titration = StrongAcidStrongBaseTitration(Ca=0.1, Va=0.05, Cb=0.1)
    Vb = np.linspace(0.0, 0.1, 50)
    result = titration.curve(Vb)
    assert result.pH.shape == result.Vb.shape == Vb.shape


def test_gran_plot_recovers_exact_equivalence_volume():
    from chemistrykit.solutions.systems.titration import gran_plot

    titration = StrongAcidStrongBaseTitration(Ca=0.080, Va=0.025, Cb=0.100)
    Ve = titration.equivalence_volume()
    Vb = np.linspace(0.3 * Ve, 0.9 * Ve, 15)
    result = gran_plot(Vb, titration.pH_at(Vb), Va=0.025)
    assert result.equivalence_volume == pytest.approx(Ve, rel=1e-6)
    assert result.slope == pytest.approx(-0.100, rel=1e-6)
    assert result.intercept == pytest.approx(0.080 * 0.025, rel=1e-6)


def test_gran_plot_is_robust_to_seeded_ph_noise():
    from chemistrykit.solutions.systems.titration import gran_plot

    rng = np.random.default_rng(12345)
    titration = StrongAcidStrongBaseTitration(Ca=0.100, Va=0.050, Cb=0.100)
    Ve = titration.equivalence_volume()
    Vb = np.linspace(0.5 * Ve, 0.95 * Ve, 20)
    noisy_pH = titration.pH_at(Vb) + rng.normal(0.0, 0.005, Vb.size)
    result = gran_plot(Vb, noisy_pH, Va=0.050)
    assert result.equivalence_volume == pytest.approx(Ve, rel=5e-3)
