"""Tests for chemistrykit.electrochem.systems.nernst against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import FARADAY, R
from chemistrykit.electrochem.systems.nernst import (
    activity_corrected_reaction_quotient,
    concentration_cell_potential,
    nernst_potential,
    nernst_potential_with_activity,
)


def test_nernst_reduces_to_standard_potential_at_unit_activity():
    """At Q=1 (unit activity), the Nernst equation gives exactly E_standard."""
    assert nernst_potential(E_standard=0.80, n=1, Q=1.0) == pytest.approx(0.80)
    assert nernst_potential(E_standard=-0.76, n=2, Q=1.0) == pytest.approx(-0.76)


def test_nernst_matches_manual_formula():
    E_standard, n, Q, T = 1.10, 2, 25.0, 310.0
    expected = E_standard - (R * T) / (n * FARADAY) * np.log(Q)
    assert nernst_potential(E_standard, n, Q, T=T) == pytest.approx(expected)


def test_nernst_potential_decreases_with_increasing_Q():
    E_standard, n = 0.34, 2
    Q = np.array([0.1, 1.0, 10.0, 100.0])
    E = nernst_potential(E_standard, n, Q)
    assert np.all(np.diff(E) < 0)


def test_concentration_cell_zero_at_equal_concentrations():
    assert concentration_cell_potential(n=1, C_cathode=0.5, C_anode=0.5) == pytest.approx(0.0, abs=1e-12)


def test_concentration_cell_positive_when_cathode_more_concentrated():
    E = concentration_cell_potential(n=1, C_cathode=0.1, C_anode=0.01)
    assert E > 0.0
    # Reversing which side is more concentrated flips the sign.
    E_rev = concentration_cell_potential(n=1, C_cathode=0.01, C_anode=0.1)
    assert E_rev == pytest.approx(-E)


def test_concentration_cell_known_59mV_per_decade():
    """A tenfold concentration ratio in a 1-electron couple at 25 degC gives ~59.2 mV."""
    E = concentration_cell_potential(n=1, C_cathode=1.0, C_anode=0.1)
    assert E == pytest.approx(0.0592, abs=2e-4)


def test_activity_corrected_quotient_reduces_to_ideal_at_zero_ionic_strength():
    c = [1e-9, 1e-9]
    z = [1, -1]
    nu = [-1.0, 1.0]
    Q_ideal = c[1] / c[0]
    Q_corrected = activity_corrected_reaction_quotient(c, z, nu)
    assert Q_corrected == pytest.approx(Q_ideal, rel=1e-3)


def test_activity_correction_deviates_from_ideal_at_higher_ionic_strength():
    """At non-negligible ionic strength, activity coefficients < 1 pull the
    activity-corrected reaction quotient away from the raw-concentration
    (ideal) value -- but only when the two species carry *different*
    charge magnitudes; equal |z| gives identical gamma for both species,
    which cancels exactly in the ratio (so that symmetric case is not a
    useful deviation test)."""
    c = [0.1, 0.1]
    z = [1, -2]
    nu = [-1.0, 1.0]
    Q_ideal = c[1] / c[0]
    Q_corrected = activity_corrected_reaction_quotient(c, z, nu)
    assert Q_corrected != pytest.approx(Q_ideal, rel=1e-6)


def test_nernst_potential_with_activity_matches_ideal_in_dilute_limit():
    c = [1e-7, 1e-13]
    z = [1, 0]
    nu = [-1.0, 1.0]
    Q_raw = c[1] / c[0]
    ideal = nernst_potential(0.50, n=1, Q=Q_raw)
    corrected = nernst_potential_with_activity(0.50, n=1, concentrations=c, charges=z, stoich_coeffs=nu)
    assert corrected == pytest.approx(ideal, abs=1e-4)
