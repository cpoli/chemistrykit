"""chemistrykit.surface: surface chemistry and heterogeneous catalysis.

Langmuir, Freundlich, and BET adsorption isotherms, each with its
standard linearization for fitting parameters from data; Langmuir-
Hinshelwood surface-reaction kinetics (single- and dual-site mechanisms);
and a turnover-number/rate-enhancement catalysis model built on
:mod:`chemistrykit.kinetics`'s Arrhenius equation.
"""

__version__ = "0.1.0"

from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.systems.bet import BETFit, BETIsotherm, bet_loading, fit_bet
from chemistrykit.surface.systems.catalysis import (
    CatalyticRateComparison,
    compare_catalyzed_rate,
    turnover_frequency,
    turnover_number,
)
from chemistrykit.surface.systems.freundlich import FreundlichFit, FreundlichIsotherm, fit_freundlich, freundlich_loading
from chemistrykit.surface.systems.langmuir import LangmuirFit, LangmuirIsotherm, fit_langmuir, langmuir_coverage
from chemistrykit.surface.systems.langmuir_hinshelwood import lh_rate_dual_site, lh_rate_single_site

__all__ = [
    "__version__",
    "AdsorptionIsotherm",
    "langmuir_coverage",
    "LangmuirIsotherm",
    "LangmuirFit",
    "fit_langmuir",
    "freundlich_loading",
    "FreundlichIsotherm",
    "FreundlichFit",
    "fit_freundlich",
    "bet_loading",
    "BETIsotherm",
    "BETFit",
    "fit_bet",
    "lh_rate_single_site",
    "lh_rate_dual_site",
    "turnover_frequency",
    "turnover_number",
    "CatalyticRateComparison",
    "compare_catalyzed_rate",
]
