"""Tests for chemistrykit.spectro.systems.atomic against the closed-form Balmer-Rydberg formula."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.spectro.systems.atomic import rydberg_wavenumber


def test_balmer_formula_equivalence():
    # Balmer (1885): lambda = B n^2 / (n^2 - 4), with B = 4 / R_H.
    n = np.arange(3, 10)
    r_h = sc.Rydberg / 100.0 / (1.0 + sc.m_e / sc.m_p)
    lam_cm = 1.0 / rydberg_wavenumber(2, n, nuclear_mass=sc.m_p)
    assert lam_cm == pytest.approx((4.0 / r_h) * n**2 / (n**2 - 4.0), rel=1e-12)


def test_h_alpha_vacuum_wavelength():
    assert 1e7 / rydberg_wavenumber(2, 3, nuclear_mass=sc.m_p) == pytest.approx(656.47, abs=0.01)


def test_nuclear_charge_scales_as_z_squared():
    assert rydberg_wavenumber(1, 2, nuclear_charge=2) == pytest.approx(4.0 * rydberg_wavenumber(1, 2))


def test_deuterium_lines_lie_above_hydrogen():
    nu_h = rydberg_wavenumber(2, 3, nuclear_mass=sc.m_p)
    nu_d = rydberg_wavenumber(2, 3, nuclear_mass=sc.physical_constants["deuteron mass"][0])
    assert nu_d > nu_h


def test_rejects_invalid_levels():
    with pytest.raises(ValueError):
        rydberg_wavenumber(3, 2)
    with pytest.raises(ValueError):
        rydberg_wavenumber(0, 2)
