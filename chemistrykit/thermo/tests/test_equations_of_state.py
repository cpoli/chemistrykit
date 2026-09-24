"""Tests for chemistrykit.thermo.systems.equations_of_state against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.thermo.systems.equations_of_state import IdealGas, RedlichKwong, VanDerWaals


def test_ideal_gas_molar_volume_at_stp():
    """The classic 22.4 L/mol at 1 atm-scale pressure and 0 degC."""
    gas = IdealGas()
    Vm = gas.molar_volume(P=101325.0, T=273.15)
    assert Vm * 1000.0 == pytest.approx(22.4, abs=0.05)


def test_ideal_gas_compressibility_factor_is_one():
    gas = IdealGas()
    assert gas.compressibility_factor(P=1.0e6, T=400.0) == pytest.approx(1.0)


def test_ideal_gas_pressure_and_molar_volume_are_inverses():
    gas = IdealGas()
    P, T = 2.5e5, 350.0
    Vm = gas.molar_volume(P, T)
    assert gas.pressure(Vm, T) == pytest.approx(P)


def test_van_der_waals_reduces_to_ideal_gas_when_a_b_zero():
    vdw = VanDerWaals(a=0.0, b=0.0)
    ideal = IdealGas()
    Vm = np.linspace(1e-3, 1e-1, 20)
    np.testing.assert_allclose(vdw.pressure(Vm, T=300.0), ideal.pressure(Vm, T=300.0))


def test_van_der_waals_molar_volume_round_trips_through_pressure():
    vdw = VanDerWaals(a=0.1448, b=3.913e-5)
    P, T = 1.0e5, 300.0
    Vm = vdw.molar_volume(P, T)
    assert vdw.pressure(Vm, T) == pytest.approx(P, rel=1e-8)


def test_van_der_waals_critical_compressibility_factor_is_three_eighths():
    """Zc = 3/8 exactly for van der Waals, independent of substance (Atkins & de Paula Table 1.6)."""
    Tc, Pc = 304.13, 7.3773e6  # CO2
    eos = VanDerWaals.from_critical_constants(Tc=Tc, Pc=Pc)
    Vc = 3.0 * eos.b
    assert Pc * Vc / (R * Tc) == pytest.approx(0.375, abs=1e-6)


def test_van_der_waals_rejects_negative_b():
    with pytest.raises(ValueError):
        VanDerWaals(a=1.0, b=-1.0)


def test_redlich_kwong_molar_volume_round_trips_through_pressure():
    rk = RedlichKwong(a=6.4239, b=2.7143e-5)
    P, T = 1.0e5, 300.0
    Vm = rk.molar_volume(P, T)
    assert rk.pressure(Vm, T) == pytest.approx(P, rel=1e-8)


def test_redlich_kwong_critical_compressibility_factor_is_one_third():
    """Zc = 1/3 exactly for Redlich-Kwong, independent of substance."""
    Tc, Pc = 304.13, 7.3773e6  # CO2
    eos = RedlichKwong.from_critical_constants(Tc=Tc, Pc=Pc)
    Vc = eos.R * Tc / (3.0 * Pc)
    assert Pc * Vc / (R * Tc) == pytest.approx(1.0 / 3.0, abs=1e-6)


def test_redlich_kwong_rejects_negative_b():
    with pytest.raises(ValueError):
        RedlichKwong(a=1.0, b=-1.0)


def test_van_der_waals_and_ideal_gas_agree_at_low_pressure():
    """At low pressure/high molar volume, real-gas corrections vanish."""
    vdw = VanDerWaals(a=0.1448, b=3.913e-5)
    ideal = IdealGas()
    P, T = 1.0e3, 500.0  # very dilute
    Vm_vdw = vdw.molar_volume(P, T)
    Vm_ideal = ideal.molar_volume(P, T)
    assert Vm_vdw == pytest.approx(Vm_ideal, rel=1e-3)


def test_van_der_waals_vapor_branch_has_larger_volume_than_liquid_branch():
    """Well below the critical point, the cubic has three real roots (liquid/unstable/vapor)."""
    Tc, Pc = 304.13, 7.3773e6
    eos = VanDerWaals.from_critical_constants(Tc=Tc, Pc=Pc)
    T, P = 0.85 * Tc, 0.5 * Pc
    vapor = eos.molar_volume(P, T, branch="vapor")
    liquid = eos.molar_volume(P, T, branch="liquid")
    assert vapor > liquid
    assert eos.pressure(vapor, T) == pytest.approx(P, rel=1e-6)
    assert eos.pressure(liquid, T) == pytest.approx(P, rel=1e-6)


def test_peng_robinson_alpha_is_one_at_critical_temperature():
    from chemistrykit.thermo.systems.equations_of_state import PengRobinson

    eos = PengRobinson(Tc=304.13, Pc=7.3773e6, omega=0.224)
    assert eos.alpha(304.13) == pytest.approx(1.0)
    assert eos.kappa == pytest.approx(0.37464 + 1.54226 * 0.224 - 0.26992 * 0.224**2)


def test_peng_robinson_critical_point_is_inflection():
    """P(Vc, Tc) = Pc and dP/dV = 0 at Vc = Zc RTc/Pc with Zc = 0.30740."""
    from chemistrykit.thermo.systems.equations_of_state import PengRobinson

    Tc, Pc = 190.56, 4.599e6
    eos = PengRobinson(Tc=Tc, Pc=Pc, omega=0.011)
    Vc = 0.307401 * R * Tc / Pc
    assert eos.pressure(Vc, Tc) == pytest.approx(Pc, rel=2e-4)
    h = 1e-3 * Vc
    slope = (eos.pressure(Vc + h, Tc) - eos.pressure(Vc - h, Tc)) / (2 * h)
    assert abs(slope * Vc / Pc) < 1e-2


def test_peng_robinson_molar_volume_round_trips_and_branches():
    from chemistrykit.thermo.systems.equations_of_state import PengRobinson

    eos = PengRobinson(Tc=304.13, Pc=7.3773e6, omega=0.224)
    P, T = 3.0e6, 270.0
    Vv = eos.molar_volume(P, T)
    Vl = eos.molar_volume(P, T, branch="liquid")
    assert Vl < Vv
    assert eos.pressure(Vv, T) == pytest.approx(P, rel=1e-8)
    assert eos.pressure(Vl, T) == pytest.approx(P, rel=1e-6)


def test_peng_robinson_approaches_ideal_gas_at_low_pressure():
    from chemistrykit.thermo.systems.equations_of_state import PengRobinson

    eos = PengRobinson(Tc=304.13, Pc=7.3773e6, omega=0.224)
    assert eos.compressibility_factor(10.0, 500.0) == pytest.approx(1.0, abs=1e-5)
