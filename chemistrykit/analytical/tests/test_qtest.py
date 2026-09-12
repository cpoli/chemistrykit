"""Tests for chemistrykit.analytical.systems.qtest: critical-value table and rejection logic."""

import pytest

from chemistrykit.analytical.systems.qtest import Q_CRITICAL_TABLE, dixon_q_test


@pytest.mark.parametrize("n", range(3, 11))
def test_q_critical_table_has_all_three_confidence_levels(n):
    assert set(Q_CRITICAL_TABLE[n]) == {0.90, 0.95, 0.99}


@pytest.mark.parametrize("n", range(3, 11))
def test_q_critical_values_increase_with_confidence(n):
    row = Q_CRITICAL_TABLE[n]
    assert row[0.90] < row[0.95] < row[0.99]


@pytest.mark.parametrize("confidence", [0.90, 0.95, 0.99])
def test_q_critical_values_decrease_with_sample_size(confidence):
    values = [Q_CRITICAL_TABLE[n][confidence] for n in range(3, 11)]
    assert all(values[i] > values[i + 1] for i in range(len(values) - 1))


def test_q_critical_table_matches_rorabacher_reference_spot_checks():
    # Rorabacher, Anal. Chem. 63, 139 (1991), Table II spot checks.
    assert Q_CRITICAL_TABLE[3][0.90] == pytest.approx(0.941)
    assert Q_CRITICAL_TABLE[3][0.99] == pytest.approx(0.994)
    assert Q_CRITICAL_TABLE[10][0.95] == pytest.approx(0.466)
    assert Q_CRITICAL_TABLE[7][0.95] == pytest.approx(0.568)


def test_dixon_q_test_rejects_clear_outlier():
    result = dixon_q_test([10.0, 10.1, 10.2, 10.3, 15.0], confidence=0.95)
    assert result.suspect_value == pytest.approx(15.0)
    assert result.reject is True


def test_dixon_q_test_retains_borderline_value():
    result = dixon_q_test([10.0, 10.1, 10.2, 10.3, 10.9], confidence=0.99)
    assert result.reject is False


def test_dixon_q_test_auto_detects_minimum_as_suspect():
    result = dixon_q_test([1.0, 10.0, 10.1, 10.2, 10.3], confidence=0.90)
    assert result.suspect_value == pytest.approx(1.0)


def test_dixon_q_test_explicit_suspect_side():
    data = [10.0, 10.1, 10.2, 10.3, 10.9]
    result_min = dixon_q_test(data, confidence=0.95, suspect="min")
    result_max = dixon_q_test(data, confidence=0.95, suspect="max")
    assert result_min.suspect_value == pytest.approx(10.0)
    assert result_max.suspect_value == pytest.approx(10.9)


def test_dixon_q_test_rejects_out_of_range_sample_size():
    with pytest.raises(ValueError):
        dixon_q_test([1.0, 2.0])
    with pytest.raises(ValueError):
        dixon_q_test(list(range(11)))


def test_dixon_q_test_rejects_invalid_confidence():
    with pytest.raises(ValueError):
        dixon_q_test([1.0, 2.0, 3.0], confidence=0.80)
