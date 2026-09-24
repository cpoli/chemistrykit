"""Tests for chemistrykit.analytical.systems.statistics: t intervals, Grubbs' test, Horwitz function."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.statistics import (
    grubbs_critical_value,
    grubbs_test,
    horrat,
    horwitz_rsd,
    t_confidence_interval,
)


@pytest.mark.parametrize(("n", "t_table"), [(2, 12.706), (3, 4.303), (5, 2.776), (11, 2.228)])
def test_t_critical_matches_tabulated_95_percent_values(n, t_table):
    ci = t_confidence_interval(np.arange(n, dtype=float), confidence=0.95)
    assert ci.t_critical == pytest.approx(t_table, abs=1e-3)


def test_t_interval_is_mean_plus_minus_t_s_over_root_n():
    data = [4.1, 4.3, 4.0, 4.4, 4.2, 4.2]
    ci = t_confidence_interval(data, confidence=0.90)
    s = np.std(data, ddof=1)
    assert ci.lower == pytest.approx(np.mean(data) - ci.t_critical * s / np.sqrt(6))
    assert ci.upper - ci.lower == pytest.approx(2 * ci.half_width)


def test_t_interval_coverage_matches_nominal_level():
    rng = np.random.default_rng(0)
    trials, hits = 4000, 0
    for _ in range(trials):
        ci = t_confidence_interval(rng.normal(5.0, 0.2, size=4), confidence=0.95)
        hits += ci.lower <= 5.0 <= ci.upper
    assert hits / trials == pytest.approx(0.95, abs=0.015)


@pytest.mark.parametrize(("n", "G_table"), [(3, 1.155), (5, 1.715), (10, 2.290), (20, 2.709)])
def test_grubbs_critical_values_match_table(n, G_table):
    assert grubbs_critical_value(n, alpha=0.05) == pytest.approx(G_table, abs=2e-3)


def test_grubbs_statistic_definition_and_rejection():
    data = [10.0, 10.1, 10.2, 10.3, 15.0]
    result = grubbs_test(data)
    assert result.G_statistic == pytest.approx(abs(15.0 - np.mean(data)) / np.std(data, ddof=1))
    assert result.reject is True
    assert grubbs_test([10.0, 10.1, 10.2, 10.3, 10.4]).reject is False


def test_grubbs_false_positive_rate_matches_alpha():
    rng = np.random.default_rng(1)
    rejections = sum(grubbs_test(rng.normal(size=8), alpha=0.05).reject for _ in range(4000))
    assert rejections / 4000 == pytest.approx(0.05, abs=0.012)


def test_horwitz_doubles_per_hundredfold_dilution():
    C = np.array([1.0, 1e-2, 1e-4, 1e-6])
    rsd = horwitz_rsd(C)
    assert rsd == pytest.approx([2.0, 4.0, 8.0, 16.0])
    assert horwitz_rsd(1e-9) == pytest.approx(2.0 * 1e-9 ** (-0.1505), rel=1e-3)


def test_horrat_is_one_at_horwitz_value():
    assert horrat(horwitz_rsd(1e-5), 1e-5) == pytest.approx(1.0)
