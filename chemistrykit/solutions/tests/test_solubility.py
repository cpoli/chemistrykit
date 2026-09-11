"""Tests for chemistrykit.solutions.systems.solubility against closed-form results."""

import pytest

from chemistrykit.solutions.systems.solubility import (
    ksp_from_molar_solubility,
    molar_solubility_from_ksp,
    molar_solubility_with_common_ion,
)


def test_ksp_and_molar_solubility_round_trip_for_1_1_salt():
    s = 1.3e-5
    Ksp = ksp_from_molar_solubility(s, cation_coeff=1, anion_coeff=1)
    assert molar_solubility_from_ksp(Ksp, cation_coeff=1, anion_coeff=1) == pytest.approx(s, rel=1e-9)


def test_ksp_and_molar_solubility_round_trip_for_1_2_salt():
    s = 2.1e-4
    Ksp = ksp_from_molar_solubility(s, cation_coeff=1, anion_coeff=2)
    assert molar_solubility_from_ksp(Ksp, cation_coeff=1, anion_coeff=2) == pytest.approx(s, rel=1e-9)


def test_ksp_formula_matches_direct_ion_product_for_caf2():
    """Ksp = [Ca2+][F-]^2 directly, cross-checked against the p^p q^q s^(p+q) formula."""
    s = 2.1e-4
    Ca_conc = s
    F_conc = 2.0 * s
    Ksp_direct = Ca_conc * F_conc**2
    Ksp_formula = ksp_from_molar_solubility(s, cation_coeff=1, anion_coeff=2)
    assert Ksp_formula == pytest.approx(Ksp_direct)


def test_common_ion_effect_suppresses_solubility():
    Ksp = 1.8e-10  # AgCl
    s_pure = molar_solubility_from_ksp(Ksp, cation_coeff=1, anion_coeff=1)
    s_common = molar_solubility_with_common_ion(Ksp, cation_coeff=1, anion_coeff=1, common_ion_conc=0.10, common_ion="anion")
    assert s_common < s_pure


def test_common_ion_effect_matches_dominant_ion_approximation():
    """When the common ion swamps the salt's own contribution, s ~= Ksp/C0."""
    Ksp = 1.8e-10
    C0 = 0.10
    s_common = molar_solubility_with_common_ion(Ksp, cation_coeff=1, anion_coeff=1, common_ion_conc=C0, common_ion="anion")
    assert s_common == pytest.approx(Ksp / C0, rel=1e-2)


def test_common_ion_effect_satisfies_exact_ksp():
    Ksp = 1.8e-10
    C0 = 0.05
    s = molar_solubility_with_common_ion(Ksp, cation_coeff=1, anion_coeff=1, common_ion_conc=C0, common_ion="cation")
    cation = s + C0
    anion = s
    assert cation * anion == pytest.approx(Ksp, rel=1e-6)


def test_common_ion_effect_symmetric_between_cation_and_anion_for_1_1_salt():
    Ksp = 2.5e-9
    C0 = 0.02
    s_cation = molar_solubility_with_common_ion(Ksp, 1, 1, C0, common_ion="cation")
    s_anion = molar_solubility_with_common_ion(Ksp, 1, 1, C0, common_ion="anion")
    assert s_cation == pytest.approx(s_anion, rel=1e-9)


def test_invalid_common_ion_argument_raises():
    with pytest.raises(ValueError):
        molar_solubility_with_common_ion(1.0e-10, 1, 1, 0.1, common_ion="neither")
