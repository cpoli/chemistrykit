"""Tests for chemistrykit.solutions.systems.acid_base against closed-form results."""

import numpy as np
import pytest

from chemistrykit.solutions.systems.acid_base import (
    Buffer,
    WeakAcid,
    WeakBase,
    h_from_ph,
    henderson_hasselbalch_ph,
    oh_from_poh,
    ph_from_h,
    poh_from_oh,
)


def test_ph_pow_are_inverses():
    assert h_from_ph(ph_from_h(3.7e-5)) == pytest.approx(3.7e-5, rel=1e-9)


def test_poh_and_oh_are_inverses():
    assert oh_from_poh(poh_from_oh(2.2e-4)) == pytest.approx(2.2e-4, rel=1e-9)


def test_neutral_water_ph_and_poh_sum_to_pkw():
    h = np.sqrt(1e-14)
    oh = 1e-14 / h
    assert ph_from_h(h) + poh_from_oh(oh) == pytest.approx(14.0)


def test_weak_acid_matches_quadratic_approximation_when_ka_ca_much_larger_than_kw():
    acid = WeakAcid(Ca=0.100, Ka=1.8e-5)
    h_exact = acid.h_concentration()
    h_approx = np.sqrt(1.8e-5 * 0.100)
    assert h_exact == pytest.approx(h_approx, rel=1e-2)


def test_weak_acid_satisfies_charge_and_mass_balance():
    acid = WeakAcid(Ca=0.050, Ka=6.3e-5)
    h = acid.h_concentration()
    oh = acid.Kw / h
    A_minus = h - oh
    HA = acid.Ca - A_minus
    # Ka = [H+][A-]/[HA]
    assert h * A_minus / HA == pytest.approx(acid.Ka, rel=1e-6)
    assert HA + A_minus == pytest.approx(acid.Ca, rel=1e-6)


def test_weak_base_is_mirror_image_of_weak_acid():
    """A WeakAcid and WeakBase with the same C and equilibrium constant give identical ion magnitudes."""
    acid = WeakAcid(Ca=0.1, Ka=1.8e-5)
    base = WeakBase(Cb=0.1, Kb=1.8e-5)
    assert acid.h_concentration() == pytest.approx(base.oh_concentration(), rel=1e-9)


def test_weak_base_ph_plus_pouh_relationship():
    base = WeakBase(Cb=0.1, Kb=1.8e-5)
    assert base.pH() + base.pOH() == pytest.approx(14.0, abs=1e-6)


def test_stronger_acid_has_greater_percent_dissociation_at_fixed_concentration():
    weak = WeakAcid(Ca=0.1, Ka=1.8e-5)
    strong = WeakAcid(Ca=0.1, Ka=1.8e-3)
    assert strong.percent_dissociation() > weak.percent_dissociation()


def test_dilution_increases_percent_dissociation():
    """Ostwald dilution law: percent dissociation increases as concentration decreases."""
    concentrated = WeakAcid(Ca=1.0, Ka=1.8e-5)
    dilute = WeakAcid(Ca=0.001, Ka=1.8e-5)
    assert dilute.percent_dissociation() > concentrated.percent_dissociation()


def test_henderson_hasselbalch_equimolar_buffer_equals_pka():
    assert henderson_hasselbalch_ph(pKa=4.76, base_conc=0.2, acid_conc=0.2) == pytest.approx(4.76)


def test_henderson_hasselbalch_matches_buffer_class():
    buf = Buffer(pKa=4.76, acid_conc=0.15, base_conc=0.05)
    expected = henderson_hasselbalch_ph(4.76, 0.05, 0.15)
    assert buf.pH() == pytest.approx(expected)


def test_buffer_from_target_ph_recovers_target():
    buf = Buffer.from_target_ph(pKa=7.21, target_pH=7.4, total_conc=0.5)
    assert buf.pH() == pytest.approx(7.4, abs=1e-9)
    assert buf.acid_conc + buf.base_conc == pytest.approx(0.5)


def test_buffer_from_target_ph_more_base_when_ph_above_pka():
    buf = Buffer.from_target_ph(pKa=7.21, target_pH=7.4, total_conc=0.5)
    assert buf.base_conc > buf.acid_conc


def test_polyprotic_fractions_sum_to_one_and_have_right_shape():
    from chemistrykit.solutions.systems.acid_base import polyprotic_fractions

    pH = np.linspace(0.0, 14.0, 57)
    alpha = polyprotic_fractions(pH, [10**-2.15, 10**-7.20, 10**-12.35])
    assert alpha.shape == (4, pH.size)
    assert np.allclose(alpha.sum(axis=0), 1.0)


def test_polyprotic_fractions_monoprotic_matches_closed_form():
    from chemistrykit.solutions.systems.acid_base import polyprotic_fractions

    Ka = 1.8e-5
    pH = np.array([3.0, 4.7447, 6.0])
    h = 10.0**-pH
    alpha = polyprotic_fractions(pH, [Ka])
    assert alpha[1] == pytest.approx(Ka / (Ka + h))
    assert alpha[0] == pytest.approx(h / (Ka + h))


def test_polyprotic_fractions_adjacent_species_equal_at_each_pka():
    """For a diprotic acid, alpha_1 = alpha_2 exactly at pH = pK2 (ratio K2/h = 1)."""
    from chemistrykit.solutions.systems.acid_base import polyprotic_fractions

    alpha = polyprotic_fractions(10.33, [10**-6.35, 10**-10.33])
    assert alpha[1, 0] == pytest.approx(alpha[2, 0], rel=1e-12)


def test_buffer_capacity_matches_numerical_derivative_of_charge_balance():
    from chemistrykit.solutions.systems.acid_base import buffer_capacity

    C, Ka, Kw = 0.05, 1.8e-5, 1.0e-14

    def base_added(pH):
        h = 10.0**-pH
        return C * Ka / (Ka + h) + Kw / h - h

    for pH in (3.5, 4.74, 6.0, 9.0, 11.5):
        dpH = 1e-6
        numeric = (base_added(pH + dpH) - base_added(pH - dpH)) / (2 * dpH)
        assert float(buffer_capacity(pH, C, Ka, Kw)) == pytest.approx(numeric, rel=1e-6)


def test_buffer_capacity_peaks_at_pka():
    from chemistrykit.solutions.systems.acid_base import buffer_capacity

    Ka = 10**-4.76
    pH = np.linspace(3.5, 6.0, 2501)
    beta = buffer_capacity(pH, 0.1, Ka)
    assert pH[np.argmax(beta)] == pytest.approx(4.76, abs=2e-3)
    assert beta.max() == pytest.approx(np.log(10) * 0.1 / 4, rel=1e-3)
