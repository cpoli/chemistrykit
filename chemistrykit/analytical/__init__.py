"""chemistrykit.analytical: analytical chemistry.

Acid-base titration curves are provided by
:mod:`chemistrykit.solutions.systems.titration`; this subpackage adds
redox and complexometric (EDTA) titration-curve simulation with
equivalence-point detection; chromatographic plate theory and the van
Deemter equation (resolution, selectivity); linear-regression calibration
curves with IUPAC-convention limits of detection/quantitation; and
propagation-of-uncertainty formulas plus Dixon's Q-test for outlier
rejection. Also: Gran plots, Kovats retention indices, Purnell's
resolution equation, Student's t confidence intervals, Grubbs' outlier
test, the Horwitz precision function, and Savitzky-Golay smoothing.
"""

__version__ = "0.1.0"

from chemistrykit.analytical.core.base_system import TitrationCurve, TitrationCurveResult
from chemistrykit.analytical.systems.calibration import LinearCalibration, fit_calibration
from chemistrykit.analytical.systems.chromatography import (
    kovats_retention_index,
    minimum_plate_height,
    optimum_flow_velocity,
    plate_height,
    purnell_resolution,
    resolution,
    retention_factor,
    selectivity_factor,
    simulate_chromatogram,
    theoretical_plates,
    van_deemter_H,
)
from chemistrykit.analytical.systems.qtest import Q_CRITICAL_TABLE, QTestResult, dixon_q_test
from chemistrykit.analytical.systems.smoothing import savitzky_golay, savitzky_golay_coefficients
from chemistrykit.analytical.systems.statistics import (
    ConfidenceIntervalResult,
    GrubbsTestResult,
    grubbs_critical_value,
    grubbs_test,
    horrat,
    horwitz_rsd,
    t_confidence_interval,
)
from chemistrykit.analytical.systems.titration import EDTATitration, GranPlotResult, RedoxTitration, gran_plot
from chemistrykit.analytical.systems.uncertainty import (
    propagate_power,
    propagate_product,
    propagate_sum,
    propagate_uncertainty,
)

__all__ = [
    "__version__",
    "TitrationCurve",
    "TitrationCurveResult",
    "RedoxTitration",
    "EDTATitration",
    "GranPlotResult",
    "gran_plot",
    "theoretical_plates",
    "plate_height",
    "van_deemter_H",
    "optimum_flow_velocity",
    "minimum_plate_height",
    "retention_factor",
    "selectivity_factor",
    "resolution",
    "simulate_chromatogram",
    "kovats_retention_index",
    "purnell_resolution",
    "LinearCalibration",
    "fit_calibration",
    "propagate_sum",
    "propagate_product",
    "propagate_power",
    "propagate_uncertainty",
    "Q_CRITICAL_TABLE",
    "QTestResult",
    "dixon_q_test",
    "ConfidenceIntervalResult",
    "t_confidence_interval",
    "GrubbsTestResult",
    "grubbs_critical_value",
    "grubbs_test",
    "horwitz_rsd",
    "horrat",
    "savitzky_golay_coefficients",
    "savitzky_golay",
]
