r"""A simplified constant-current battery discharge/capacity model.

**Simplified model, flagged explicitly**: see
:class:`chemistrykit.electrochem.core.base_system.BatteryDischargeModel`
for what is and is not captured here. See Linden & Reddy, *Handbook of
Batteries*, 3rd ed., Ch. 3.3, for Peukert's law, and the original: W.
Peukert, *Elektrotechnische Zeitschrift* 20, 20 (1897).
"""

from __future__ import annotations

import numpy as np

from chemistrykit.electrochem.core.base_system import BatteryDischargeModel

__all__ = ["peukert_discharge_time", "effective_capacity", "ConstantCurrentBattery"]


def peukert_discharge_time(capacity_peukert: float, current: float, k: float = 1.0) -> float:
    r"""Peukert's-law discharge time at constant current: :math:`t = C_p / I^k`.

    :math:`C_p` (the "Peukert capacity", in :math:`A^k\cdot h`) and the
    Peukert exponent :math:`k` are empirical constants of a given
    battery chemistry/construction, fitted from discharge tests at
    several currents (Linden & Reddy, *Handbook of Batteries*, 3rd ed.,
    Ch. 3.3, eq. 3.5). :math:`k=1` is the ideal case (discharge time
    inversely proportional to current, i.e. a current-independent
    delivered capacity, see :func:`effective_capacity`); real batteries
    typically have :math:`k>1` (lead-acid: :math:`k\approx1.1-1.3`),
    reflecting the fact that a higher discharge rate delivers *less*
    total capacity than a lower one before the same cutoff voltage is
    reached.

    Parameters
    ----------
    capacity_peukert : float
        Peukert capacity constant :math:`C_p`, in :math:`A^k\cdot h`
        (numerically equal to the rated capacity in Ah only when `k=1`).
    current : float or array-like of float
        Constant discharge current, in A.
    k : float, default 1.0
        Peukert exponent, :math:`\geq1`.

    Returns
    -------
    float or ndarray
        Discharge time, in h.

    Examples
    --------
    At `k=1` (the ideal case), discharge time is exactly inversely
    proportional to current:

    >>> round(float(peukert_discharge_time(capacity_peukert=10.0, current=2.0, k=1.0)), 4)
    5.0
    >>> round(float(peukert_discharge_time(capacity_peukert=10.0, current=4.0, k=1.0)), 4)
    2.5

    At `k>1`, doubling the current more than halves the discharge time:

    >>> t1 = peukert_discharge_time(10.0, current=2.0, k=1.2)
    >>> t2 = peukert_discharge_time(10.0, current=4.0, k=1.2)
    >>> bool(t2 < t1 / 2.0)
    True
    """
    current = np.asarray(current, dtype=np.float64)
    result = capacity_peukert / current**k
    return float(result) if result.ndim == 0 else result


def effective_capacity(capacity_peukert: float, current: float, k: float = 1.0):
    r"""Effective delivered capacity at a given discharge current: :math:`C_{eff}(I) = I\,t(I) = C_p I^{1-k}`.

    At the ideal Peukert exponent :math:`k=1`, this is exactly
    :math:`C_p` regardless of current -- a battery that delivers its
    full rated capacity no matter how fast it is discharged. Real
    batteries (:math:`k>1`) deliver progressively less effective
    capacity as the discharge current increases, the practical meaning
    of "capacity fade with rate" (Linden & Reddy, *Handbook of
    Batteries*, 3rd ed., Ch. 3.3).

    Parameters
    ----------
    capacity_peukert : float
        Peukert capacity constant, in :math:`A^k\cdot h`.
    current : float or array-like of float
        Constant discharge current, in A.
    k : float, default 1.0
        Peukert exponent.

    Returns
    -------
    float or ndarray
        Effective capacity, in Ah.

    Examples
    --------
    At `k=1`, effective capacity is independent of current (the ideal,
    rate-independent battery):

    >>> c_low = effective_capacity(10.0, current=1.0, k=1.0)
    >>> c_high = effective_capacity(10.0, current=5.0, k=1.0)
    >>> round(c_low, 6), round(c_high, 6)
    (10.0, 10.0)

    At `k>1`, effective capacity decreases as discharge current
    increases -- the real-battery "capacity fade at high rate" effect:

    >>> c_low = effective_capacity(10.0, current=1.0, k=1.2)
    >>> c_high = effective_capacity(10.0, current=5.0, k=1.2)
    >>> bool(c_high < c_low)
    True
    """
    current = np.asarray(current, dtype=np.float64)
    result = np.asarray(current * peukert_discharge_time(capacity_peukert, current, k), dtype=np.float64)
    return float(result) if result.ndim == 0 else result


class ConstantCurrentBattery(BatteryDischargeModel):
    r"""A battery discharged at constant current, with a Peukert-law-limited runtime and an ohmic voltage sag.

    **Simplified model** (see the module and
    :class:`~chemistrykit.electrochem.core.base_system.BatteryDischargeModel`
    docstrings): the open-circuit voltage is treated as constant at
    `v_nominal` until the Peukert-law discharge time is reached, at
    which point the battery is treated as instantaneously exhausted
    (voltage drops to 0); the only voltage variation modeled during
    discharge is a constant ohmic sag :math:`IR_{internal}`, not the
    smoothly declining open-circuit-voltage-vs-state-of-charge curve of
    a real cell.

    Parameters
    ----------
    capacity_peukert : float
        Peukert capacity constant, in :math:`A^k\cdot h`.
    current : float
        Constant discharge current, in A.
    v_nominal : float
        Open-circuit (fully charged) terminal voltage, in V.
    internal_resistance : float, default 0.0
        Internal (ohmic) resistance, in ohms.
    k : float, default 1.0
        Peukert exponent.

    Examples
    --------
    >>> battery = ConstantCurrentBattery(capacity_peukert=2.0, current=1.0, v_nominal=3.7, internal_resistance=0.05)
    >>> round(battery.discharge_time(), 4)
    2.0
    >>> round(float(battery.terminal_voltage(1.0)), 4)
    3.65
    >>> round(float(battery.terminal_voltage(3.0)), 4)
    0.0
    """

    def __init__(self, capacity_peukert: float, current: float, v_nominal: float, internal_resistance: float = 0.0, k: float = 1.0):
        self.capacity_peukert = capacity_peukert
        self.current = current
        self.v_nominal = v_nominal
        self.internal_resistance = internal_resistance
        self.k = k

    def discharge_time(self) -> float:
        """Total discharge time (runtime until exhausted), in h.

        Returns
        -------
        float
        """
        return peukert_discharge_time(self.capacity_peukert, self.current, self.k)

    def state_of_charge(self, t):
        t = np.asarray(t, dtype=np.float64)
        t_total = self.discharge_time()
        result = np.clip(1.0 - t / t_total, 0.0, 1.0)
        return float(result) if result.ndim == 0 else result

    def terminal_voltage(self, t):
        t = np.asarray(t, dtype=np.float64)
        soc = self.state_of_charge(t)
        v_under_load = self.v_nominal - self.current * self.internal_resistance
        result = np.where(soc > 0.0, v_under_load, 0.0)
        return float(result) if result.ndim == 0 else result
