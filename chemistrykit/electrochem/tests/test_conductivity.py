"""Tests for chemistrykit.electrochem.systems.conductivity against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.electrochem.systems.conductivity import (
    LIMITING_IONIC_CONDUCTIVITIES,
    fit_kohlrausch_law,
    kohlrausch_molar_conductivity,
    limiting_molar_conductivity,
)


def test_limiting_molar_conductivity_is_stoichiometric_sum():
    T = LIMITING_IONIC_CONDUCTIVITIES
    assert limiting_molar_conductivity({"Mg2+": 1, "Cl-": 2}) == pytest.approx(T["Mg2+"] + 2 * T["Cl-"])


def test_acetic_acid_from_strong_electrolytes():
    """Kohlrausch's classic trick: L(HAc) = L(HCl) + L(NaAc) - L(NaCl)."""
    hcl = limiting_molar_conductivity({"H+": 1, "Cl-": 1})
    naac = limiting_molar_conductivity({"Na+": 1, "CH3COO-": 1})
    nacl = limiting_molar_conductivity({"Na+": 1, "Cl-": 1})
    direct = limiting_molar_conductivity({"H+": 1, "CH3COO-": 1})
    assert hcl + naac - nacl == pytest.approx(direct)
    assert direct * 1e4 == pytest.approx(390.5, abs=0.1)


def test_square_root_law_values():
    assert kohlrausch_molar_conductivity(0.0, 0.0126, 0.0089) == pytest.approx(0.0126)
    assert kohlrausch_molar_conductivity(0.04, 0.0126, 0.0089) == pytest.approx(0.0126 - 0.0089 * 0.2)


def test_fit_recovers_parameters_from_noisy_data():
    rng = np.random.default_rng(0)
    c = np.linspace(1e-4, 1e-2, 20)
    data = kohlrausch_molar_conductivity(c, 0.0150, 0.0095) + rng.normal(0.0, 1e-6, c.size)
    fit = fit_kohlrausch_law(c, data)
    assert fit.limiting_molar_conductivity == pytest.approx(0.0150, rel=1e-3)
    assert fit.kohlrausch_coefficient == pytest.approx(0.0095, rel=2e-2)
    assert fit.r_squared > 0.999
