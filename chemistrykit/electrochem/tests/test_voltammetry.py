"""Tests for chemistrykit.electrochem.systems.voltammetry against closed-form/known results."""

import numpy as np
import pytest
from scipy.integrate import quad

from chemistrykit.constants import FARADAY, R
from chemistrykit.electrochem.systems.voltammetry import (
    cottrell_current,
    ilkovic_diffusion_current,
    polarographic_wave_current,
    randles_sevcik_peak_current,
)


def test_cottrell_charge_matches_anson_integral():
    """Integrating the Cottrell current gives Q = 2nFAC sqrt(Dt/pi)."""
    n, A, C, D, t = 1, 1e-6, 2.0, 1e-9, 5.0
    Q, _ = quad(lambda s: cottrell_current(s, n, A, C, D), 0.0, t)
    assert Q == pytest.approx(2 * n * FARADAY * A * C * np.sqrt(D * t / np.pi), rel=1e-6)


def test_ilkovic_scaling_exponents():
    base = ilkovic_diffusion_current(1, 1e-5, 2.0, 4.0, 1.0)
    assert ilkovic_diffusion_current(1, 4e-5, 2.0, 4.0, 1.0) == pytest.approx(2.0 * base)
    assert ilkovic_diffusion_current(1, 1e-5, 16.0, 4.0, 1.0) == pytest.approx(4.0 * base)
    assert ilkovic_diffusion_current(1, 1e-5, 2.0, 256.0, 1.0) == pytest.approx(2.0 * base)
    assert ilkovic_diffusion_current(1, 1e-5, 2.0, 4.0, 1.0, average=True) == pytest.approx(base * 607 / 708)


def test_polarographic_wave_shape():
    E_half, n, i_d = -0.5, 2, 10.0
    assert polarographic_wave_current(E_half, E_half, i_d, n) == pytest.approx(i_d / 2)
    assert polarographic_wave_current(-1.5, E_half, i_d, n) == pytest.approx(i_d)
    assert polarographic_wave_current(0.5, E_half, i_d, n) == pytest.approx(0.0, abs=1e-10)
    E = np.array([-0.55, -0.45])
    i = polarographic_wave_current(E, E_half, i_d, n)
    recovered = E_half + R * 298.15 / (n * FARADAY) * np.log((i_d - i) / i)
    assert recovered == pytest.approx(E)


def test_randles_sevcik_matches_textbook_constant():
    """At 25 degC, i_p = 2.69e5 n^1.5 A D^0.5 C v^0.5 in cgs units (A cm^2, D cm^2/s, C mol/cm^3)."""
    n, A_cm2, D_cm2, C_molcm3, v = 2, 0.1, 1e-5, 1e-6, 0.05
    ip = randles_sevcik_peak_current(v, n, A_cm2 * 1e-4, C_molcm3 * 1e6, D_cm2 * 1e-4)
    assert ip == pytest.approx(2.69e5 * n**1.5 * A_cm2 * np.sqrt(D_cm2) * C_molcm3 * np.sqrt(v), rel=2e-3)
