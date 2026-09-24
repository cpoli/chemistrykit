"""Tests for chemistrykit.photochem.systems.electron_transfer."""

import numpy as np
import pytest

from chemistrykit.constants import ELEMENTARY_CHARGE, K_B
from chemistrykit.photochem.systems.electron_transfer import rehm_weller_free_energy, rehm_weller_quenching_rate


def test_free_energy():
    assert rehm_weller_free_energy(1.2, -1.9, 3.3, work_term=-0.06) == pytest.approx(-0.26)


def test_quenching_rate_limits_and_monotonic():
    kT = K_B * 298.15 / ELEMENTARY_CHARGE
    assert rehm_weller_quenching_rate(0.0) == pytest.approx(2e10 / (1 + 0.25 * (np.exp(0.104 / kT) + 1)))
    assert rehm_weller_quenching_rate(-100.0) == pytest.approx(1.6e10, rel=1e-3)
    dG = np.linspace(-2.0, 1.0, 50)
    kq = rehm_weller_quenching_rate(dG)
    assert np.all(np.diff(kq) < 0)  # no Marcus inverted region
    assert kq[-1] < 1e-6 * kq[0]
