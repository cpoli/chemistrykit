"""Abstract base class for battery discharge models, and the result container.

Most of :mod:`chemistrykit.electrochem` -- the Nernst equation, standard
reduction potentials and cell-potential combination, Butler-Volmer
kinetics, and Faraday's laws of electrolysis -- is a set of closed-form
algebraic relationships with no shared polymorphic behavior worth
capturing in an ABC, so (following
``chemistrykit.thermo.core.base_system``'s precedent for its reaction/
phase-equilibrium functions, and ``chemistrykit.kinetics.systems.arrhenius``/
``enzyme``'s before that) those live directly as functions and small
result dataclasses in their own ``systems/`` modules.

:mod:`chemistrykit.electrochem.systems.battery`'s discharge/capacity
models are the one genuine exception: a battery's terminal voltage and
remaining state of charge are both functions of elapsed discharge time,
:math:`t`, and more than one such model exists (an ideal constant-current
model and a Peukert's-law rate-dependent one), so
:class:`BatteryDischargeModel` captures the one shared interface, the
same way :class:`chemistrykit.thermo.core.base_system.EquationOfState`
does for the three pure-substance equations of state.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

__all__ = ["DischargeResult", "BatteryDischargeModel"]


@dataclass
class DischargeResult:
    """Container for the output of a :meth:`BatteryDischargeModel.discharge_curve` call.

    Mirrors :class:`chemistrykit.solutions.core.base_system.TitrationResult`
    (a stable, dataclass return type for a curve traced out over an
    independent variable -- there, titrant volume; here, elapsed time).
    """

    t: np.ndarray
    """ndarray: Elapsed discharge time, in the same time unit as the model's rate constants (h if `current` is in A and `capacity` in Ah)."""

    voltage: np.ndarray
    """ndarray: Terminal voltage at each time in `t`, in V."""

    state_of_charge: np.ndarray
    """ndarray: Fraction of rated capacity remaining at each time in `t`, in [0, 1]."""


class BatteryDischargeModel(ABC):
    """Common interface for a simplified constant-current battery discharge model.

    **Simplified model, flagged explicitly**: a real battery's discharge
    behavior involves a nonlinear open-circuit-voltage-vs-state-of-charge
    curve, temperature dependence, cycle aging, and relaxation effects
    under load -- none of that is modeled here. This captures only the
    two textbook-level effects usually taught alongside each other: a
    constant ohmic (IR) voltage drop under load, and the empirical
    Peukert's-law dependence of *effective* capacity on discharge rate
    (W. Peukert, *Elektrotechnische Zeitschrift* 20, 20 (1897); see
    Linden & Reddy, *Handbook of Batteries*, 3rd ed., Ch. 3.3).

    Concrete subclasses implement :meth:`state_of_charge` and
    :meth:`terminal_voltage`; :meth:`discharge_curve` is then available on
    every subclass for free, mirroring how
    :meth:`chemistrykit.solutions.core.base_system.Titration.curve` is
    shared machinery built once atop each subclass's ``pH_at``.
    """

    @abstractmethod
    def state_of_charge(self, t):
        """Return the fraction of rated capacity remaining at time(s) `t`.

        Parameters
        ----------
        t : float or array-like of float
            Elapsed discharge time.

        Returns
        -------
        float or ndarray
            In :math:`[0, 1]`.
        """

    @abstractmethod
    def terminal_voltage(self, t):
        """Return the terminal voltage at time(s) `t`.

        Parameters
        ----------
        t : float or array-like of float
            Elapsed discharge time.

        Returns
        -------
        float or ndarray
            Terminal voltage, in V (0 once the battery is exhausted).
        """

    def discharge_curve(self, t) -> DischargeResult:
        """Compute the full discharge curve over a range of elapsed times.

        Parameters
        ----------
        t : array-like of float
            Elapsed discharge times to evaluate at.

        Returns
        -------
        DischargeResult
        """
        t = np.atleast_1d(np.asarray(t, dtype=np.float64))
        voltage = np.asarray(self.terminal_voltage(t), dtype=np.float64)
        soc = np.asarray(self.state_of_charge(t), dtype=np.float64)
        return DischargeResult(t=t, voltage=voltage, state_of_charge=soc)
