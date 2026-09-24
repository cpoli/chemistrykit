r"""Electrolytic conductivity: Kohlrausch's law of independent migration and square-root law.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Topic 16C
(molar conductivity, Kohlrausch's laws, and the limiting ionic
conductivities tabulated below), or Bard & Faulkner, *Electrochemical
Methods*, 2nd ed., Ch. 2.3. The original compilation: F. Kohlrausch and
L. Holborn, *Das Leitvermögen der Elektrolyte* (Teubner, Leipzig, 1898).

All conductivities are in SI units: molar conductivity in
S m\ :sup:`2` mol\ :sup:`-1`, concentration in mol/L (the square-root law
is conventionally written against molarity).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.electrochem.utils.regression import LinearFit, linear_fit

__all__ = [
    "LIMITING_IONIC_CONDUCTIVITIES",
    "limiting_molar_conductivity",
    "kohlrausch_molar_conductivity",
    "KohlrauschFit",
    "fit_kohlrausch_law",
]

#: Limiting (infinite-dilution) molar ionic conductivities
#: :math:`\lambda_i^\circ` in water at 298.15 K, in S m^2 mol^-1, per mole
#: of the ion as written (Atkins & de Paula, *Physical Chemistry*, 11th
#: ed., Table 16C.2).
LIMITING_IONIC_CONDUCTIVITIES: dict[str, float] = {
    "H+": 34.96e-3,
    "Li+": 3.87e-3,
    "Na+": 5.01e-3,
    "K+": 7.35e-3,
    "Ag+": 6.19e-3,
    "Mg2+": 10.60e-3,
    "Ca2+": 11.90e-3,
    "OH-": 19.91e-3,
    "Cl-": 7.63e-3,
    "Br-": 7.81e-3,
    "NO3-": 7.14e-3,
    "CH3COO-": 4.09e-3,
    "SO4^2-": 16.00e-3,
}


def limiting_molar_conductivity(ions: dict[str, float], table: dict[str, float] | None = None) -> float:
    r"""Kohlrausch's law of independent migration: :math:`\Lambda_m^\circ = \sum_i \nu_i\lambda_i^\circ`.

    At infinite dilution every ion migrates independently of its
    counter-ions, so an electrolyte's limiting molar conductivity is the
    stoichiometry-weighted sum of per-ion contributions (Atkins & de
    Paula, *Physical Chemistry*, 11th ed., Topic 16C.1(b)).

    Parameters
    ----------
    ions : dict of str to float
        Map from ion name (a key of `table`) to the number of those ions
        per formula unit, e.g. ``{"Mg2+": 1, "Cl-": 2}`` for MgCl2.
        Negative multiples are allowed, which is how Kohlrausch's law is
        used to *combine* measured electrolytes.
    table : dict of str to float, optional
        Limiting ionic conductivities, in S m^2 mol^-1; defaults to
        :data:`LIMITING_IONIC_CONDUCTIVITIES`.

    Returns
    -------
    float
        Limiting molar conductivity, in S m^2 mol^-1.

    Examples
    --------
    >>> round(limiting_molar_conductivity({"Na+": 1, "Cl-": 1}) * 1e4, 1)  # S cm^2/mol
    126.4

    Independent migration means the *difference* between two salts with a
    common anion depends only on the cations:

    >>> kcl = limiting_molar_conductivity({"K+": 1, "Cl-": 1})
    >>> kno3 = limiting_molar_conductivity({"K+": 1, "NO3-": 1})
    >>> nacl = limiting_molar_conductivity({"Na+": 1, "Cl-": 1})
    >>> nano3 = limiting_molar_conductivity({"Na+": 1, "NO3-": 1})
    >>> round(kcl - nacl, 8) == round(kno3 - nano3, 8)
    True
    """
    table = LIMITING_IONIC_CONDUCTIVITIES if table is None else table
    return float(sum(nu * table[name] for name, nu in ions.items()))


def kohlrausch_molar_conductivity(c, Lambda0: float, K: float):
    r"""Kohlrausch's square-root law: :math:`\Lambda_m = \Lambda_m^\circ - K\sqrt{c}`.

    Empirically valid for strong electrolytes at low concentration
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Topic 16C.1(b));
    Debye, Hückel and Onsager later derived the coefficient `K` from
    ion-atmosphere theory.

    Parameters
    ----------
    c : float or array-like of float
        Molar concentration, in mol/L.
    Lambda0 : float
        Limiting molar conductivity, in S m^2 mol^-1.
    K : float
        Kohlrausch coefficient, in S m^2 mol^-1 (mol/L)^-1/2.

    Returns
    -------
    float or ndarray
        Molar conductivity, in S m^2 mol^-1.

    Examples
    --------
    >>> round(kohlrausch_molar_conductivity(0.0, Lambda0=0.0126, K=0.0089), 6)
    0.0126
    >>> round(float(kohlrausch_molar_conductivity(0.01, Lambda0=0.0126, K=0.0089)), 6)
    0.01171
    """
    result = Lambda0 - K * np.sqrt(np.asarray(c, dtype=np.float64))
    return float(result) if result.ndim == 0 else result


@dataclass
class KohlrauschFit:
    r"""Result of fitting :math:`\Lambda_m` vs. :math:`\sqrt{c}` data to Kohlrausch's square-root law."""

    limiting_molar_conductivity: float
    r"""float: Fitted :math:`\Lambda_m^\circ` (the intercept at :math:`c\to0`), in S m^2 mol^-1."""

    kohlrausch_coefficient: float
    """float: Fitted coefficient `K` (minus the slope)."""

    r_squared: float
    """float: Coefficient of determination of the linear fit."""


def fit_kohlrausch_law(c, Lambda_m) -> KohlrauschFit:
    r"""Extrapolate measured molar conductivities to infinite dilution.

    Fits :math:`\Lambda_m` against :math:`\sqrt{c}` by ordinary least
    squares; the intercept is :math:`\Lambda_m^\circ` and minus the slope
    is `K` -- exactly Kohlrausch's graphical extrapolation.

    Parameters
    ----------
    c : array-like of float
        Molar concentrations, in mol/L.
    Lambda_m : array-like of float
        Measured molar conductivities, in S m^2 mol^-1.

    Returns
    -------
    KohlrauschFit

    Examples
    --------
    >>> import numpy as np
    >>> c = np.array([1e-4, 1e-3, 5e-3, 1e-2])
    >>> fit = fit_kohlrausch_law(c, kohlrausch_molar_conductivity(c, 0.0126, 0.0089))
    >>> round(fit.limiting_molar_conductivity, 6), round(fit.kohlrausch_coefficient, 6)
    (0.0126, 0.0089)
    """
    fit: LinearFit = linear_fit(np.sqrt(np.asarray(c, dtype=np.float64)), np.asarray(Lambda_m, dtype=np.float64))
    return KohlrauschFit(limiting_molar_conductivity=float(fit.intercept), kohlrausch_coefficient=float(-fit.slope), r_squared=fit.r_squared)
