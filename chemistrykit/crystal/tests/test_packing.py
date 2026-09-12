"""Tests for chemistrykit.crystal.systems.packing against exact closed-form packing fractions."""

import numpy as np
import pytest

from chemistrykit.crystal.systems.packing import (
    BodyCenteredCubicPacking,
    FaceCenteredCubicPacking,
    HexagonalClosePacking,
    SimpleCubicPacking,
)


def test_simple_cubic_packing_fraction_is_pi_over_6():
    assert SimpleCubicPacking().packing_fraction() == pytest.approx(np.pi / 6.0)


def test_simple_cubic_coordination_and_atoms_per_cell():
    sc = SimpleCubicPacking()
    assert sc.coordination_number == 6
    assert sc.atoms_per_cell == 1.0


def test_bcc_packing_fraction_is_sqrt3_pi_over_8():
    assert BodyCenteredCubicPacking().packing_fraction() == pytest.approx(np.sqrt(3.0) * np.pi / 8.0)


def test_bcc_coordination_and_atoms_per_cell():
    bcc = BodyCenteredCubicPacking()
    assert bcc.coordination_number == 8
    assert bcc.atoms_per_cell == 2.0


def test_fcc_packing_fraction_is_pi_over_3sqrt2():
    assert FaceCenteredCubicPacking().packing_fraction() == pytest.approx(np.pi / (3.0 * np.sqrt(2.0)))


def test_fcc_coordination_and_atoms_per_cell():
    fcc = FaceCenteredCubicPacking()
    assert fcc.coordination_number == 12
    assert fcc.atoms_per_cell == 4.0


def test_hcp_ideal_packing_fraction_matches_fcc():
    hcp = HexagonalClosePacking()
    fcc = FaceCenteredCubicPacking()
    assert hcp.packing_fraction() == pytest.approx(fcc.packing_fraction(), rel=1e-9)


def test_hcp_coordination_number_matches_fcc():
    assert HexagonalClosePacking().coordination_number == FaceCenteredCubicPacking().coordination_number == 12


def test_hcp_ideal_c_over_a_is_sqrt_8_over_3():
    assert HexagonalClosePacking().c_over_a == pytest.approx(np.sqrt(8.0 / 3.0))


def test_hcp_nonideal_c_over_a_reduces_packing_fraction():
    ideal = HexagonalClosePacking().packing_fraction()
    zinc_like = HexagonalClosePacking(c_over_a=1.856).packing_fraction()
    assert zinc_like < ideal


def test_packing_fraction_independent_of_lattice_constant():
    fcc = FaceCenteredCubicPacking()
    assert fcc.packing_fraction(a=1.0) == pytest.approx(fcc.packing_fraction(a=5.0))


def test_packing_fraction_ordering_sc_lt_bcc_lt_fcc():
    sc = SimpleCubicPacking().packing_fraction()
    bcc = BodyCenteredCubicPacking().packing_fraction()
    fcc = FaceCenteredCubicPacking().packing_fraction()
    assert sc < bcc < fcc
