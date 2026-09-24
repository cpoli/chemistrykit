"""Tests for chemistrykit.thermo.systems.fugacity against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equations_of_state import IdealGas, PengRobinson, RedlichKwong, VanDerWaals
from chemistrykit.thermo.systems.fugacity import fugacity, fugacity_coefficient, saturation_pressure

TC, PC = 304.13, 7.3773e6


def test_ideal_gas_fugacity_equals_pressure():
    assert fugacity_coefficient(IdealGas(), 3.0e6, 350.0) == pytest.approx(1.0, abs=1e-12)
    assert fugacity(IdealGas(), 3.0e6, 350.0) == pytest.approx(3.0e6)


@pytest.mark.parametrize("branch,P,T", [("vapor", 1.0e6, 300.0), ("vapor", 2.0e7, 400.0), ("liquid", 5.0e6, 260.0)])
def test_van_der_waals_fugacity_matches_closed_form(branch, P, T):
    vdw = VanDerWaals.from_critical_constants(TC, PC)
    V = vdw.molar_volume(P, T, branch=branch)
    RT = R * T
    ln_phi = vdw.b / (V - vdw.b) - 2 * vdw.a / (RT * V) - np.log(P * (V - vdw.b) / RT)
    assert np.log(fugacity_coefficient(vdw, P, T, branch)) == pytest.approx(ln_phi, rel=1e-9, abs=1e-12)


def test_redlich_kwong_fugacity_matches_closed_form():
    rk = RedlichKwong.from_critical_constants(TC, PC)
    P, T = 4.0e6, 350.0
    V = rk.molar_volume(P, T)
    Z = P * V / (R * T)
    A = rk.a * P / (R**2 * T**2.5)
    B = rk.b * P / (R * T)
    ln_phi = Z - 1 - np.log(Z - B) - A / B * np.log(1 + B / Z)
    assert np.log(fugacity_coefficient(rk, P, T)) == pytest.approx(ln_phi, rel=1e-9)


def test_low_pressure_limit_is_second_virial():
    """ln(phi) -> B2 P/(RT) with B2 = b - a/(RT) for van der Waals."""
    vdw = VanDerWaals.from_critical_constants(TC, PC)
    P, T = 1.0e3, 400.0
    B2 = vdw.b - vdw.a / (R * T)
    assert np.log(fugacity_coefficient(vdw, P, T)) == pytest.approx(B2 * P / (R * T), rel=1e-3)


def test_van_der_waals_reduced_coexistence_pressure():
    """Tabulated vdW coexistence: Tr = 0.9 -> Pr = 0.6470 (substance independent)."""
    for Tc, Pc in [(TC, PC), (150.7, 4.863e6)]:
        vdw = VanDerWaals.from_critical_constants(Tc, Pc)
        assert saturation_pressure(vdw, 0.9 * Tc) / Pc == pytest.approx(0.6470, abs=5e-4)


def test_saturation_equalizes_liquid_and_vapor_fugacities():
    pr = PengRobinson(TC, PC, 0.224)
    Ps = saturation_pressure(pr, 270.0)
    assert fugacity(pr, Ps, 270.0, "liquid") == pytest.approx(fugacity(pr, Ps, 270.0, "vapor"), rel=1e-8)


def test_peng_robinson_reproduces_acentric_factor():
    """Pitzer's definition: omega = -1 - log10(Psat/Pc) at Tr = 0.7."""
    omega = 0.224
    pr = PengRobinson(TC, PC, omega)
    Ps = saturation_pressure(pr, 0.7 * TC)
    assert -1.0 - np.log10(Ps / PC) == pytest.approx(omega, abs=0.01)


def test_saturation_pressure_rejects_supercritical_temperature():
    vdw = VanDerWaals.from_critical_constants(TC, PC)
    with pytest.raises(ValueError):
        saturation_pressure(vdw, 1.1 * TC)
