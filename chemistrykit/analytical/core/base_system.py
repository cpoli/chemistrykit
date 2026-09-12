r"""Abstract base class for potentiometric-style titration-curve models, and a note on scope.

:class:`TitrationCurve` captures the one genuinely polymorphic shape in
this domain: several titration *types* (redox, complexometric --
:mod:`chemistrykit.analytical.systems.titration`) that all reduce to
"compute some scalar response as a function of titrant volume, then find
the equivalence point as the point of steepest response change" -- the
same shape that :class:`chemistrykit.solutions.core.base_system.Titration`
already captures for acid-base titrations (pH vs. volume).

**Why a new ABC here, rather than reusing** :class:`chemistrykit.solutions.core.base_system.Titration`:
that class's abstract method is named and documented specifically as
``pH_at`` (a *hydrogen-ion* activity readout), which is the wrong
abstraction for a redox titration's electrode potential `E` or a
complexometric titration's `pM`. Rather than force those into a
pH-flavored interface (or reach across a domain boundary the way
``chemistrykit.photochem`` reaches into ``chemistrykit.kinetics`` for
substantial shared machinery -- not warranted here, since the shared
logic is a handful of lines), :class:`TitrationCurve` reimplements the
same small pattern (curve, steepest-ascent-or-descent equivalence-point
detection) with a response-type-neutral abstract method name, exactly
the way :mod:`chemistrykit.surface.utils.regression`,
:mod:`chemistrykit.electrochem.utils.regression`, and
:mod:`chemistrykit.photochem.utils.regression` each keep their own small
``linear_fit`` rather than importing one another's. Acid-base titration
curves are *not* reimplemented here at all -- :mod:`chemistrykit.analytical`'s
examples and tests use
:mod:`chemistrykit.solutions.systems.titration`'s classes directly for
that case, side by side with the new redox/complexometric models below.

Chromatography (:mod:`chemistrykit.analytical.systems.chromatography`),
calibration curves (:mod:`chemistrykit.analytical.systems.calibration`),
uncertainty propagation (:mod:`chemistrykit.analytical.systems.uncertainty`),
and the Q-test (:mod:`chemistrykit.analytical.systems.qtest`) are each a
self-contained set of formulas with no swappable sibling, so -- following
``chemistrykit.electrochem``/``chemistrykit.photochem``/``chemistrykit.surface``'s
precedent -- they stay as plain functions (plus small result dataclasses)
in their own ``systems/`` modules.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

__all__ = ["TitrationCurveResult", "TitrationCurve"]


@dataclass
class TitrationCurveResult:
    """Container for the output of a :meth:`TitrationCurve.curve` call."""

    V: np.ndarray
    """ndarray: Volume(s) of titrant added, in L (or any consistent volume unit)."""

    response: np.ndarray
    """ndarray: The titration's response variable at each volume in `V`
    (e.g. electrode potential `E` in volts, or `pM`)."""


class TitrationCurve(ABC):
    """Common base for a potentiometric titration-curve model whose response is not pH.

    Concrete subclasses implement :meth:`response_at`; :meth:`curve` and
    :meth:`find_equivalence_point` are then available for free, mirroring
    :meth:`chemistrykit.solutions.core.base_system.Titration.curve`/
    :meth:`~chemistrykit.solutions.core.base_system.Titration.find_equivalence_point`.
    """

    @abstractmethod
    def response_at(self, V: np.ndarray) -> np.ndarray:
        """Return the titration's response variable at each titrant volume in `V`.

        Parameters
        ----------
        V : ndarray
            Volume(s) of titrant added, in L.

        Returns
        -------
        ndarray
        """

    def curve(self, V) -> TitrationCurveResult:
        """Compute the full titration curve over a range of titrant volumes.

        Parameters
        ----------
        V : array-like of float
            Volumes of titrant added, in L.

        Returns
        -------
        TitrationCurveResult
        """
        V = np.atleast_1d(np.asarray(V, dtype=np.float64))
        response = np.asarray(self.response_at(V), dtype=np.float64)
        return TitrationCurveResult(V=V, response=response)

    def find_equivalence_point(self, V) -> float:
        r"""Numerically locate the equivalence point as the point of steepest response change.

        Parameters
        ----------
        V : array-like of float
            A sufficiently fine grid of titrant volumes spanning the
            equivalence point, in L.

        Returns
        -------
        float
        """
        result = self.curve(V)
        d_response_dV = np.gradient(result.response, result.V)
        idx = int(np.argmax(np.abs(d_response_dV)))
        return float(result.V[idx])
