r"""Ionic strength and Debye-Huckel activity-coefficient laws.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 5.9, or Harris,
*Quantitative Chemical Analysis*, 9th ed., Ch. 8.2; the original: P.
Debye & E. Huckel, *Phys. Z.* 24, 185 (1923).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "DEBYE_HUCKEL_A_25C",
    "ionic_strength",
    "activity_coefficient_debye_huckel_limiting",
    "activity_coefficient_debye_huckel_extended",
    "activity_coefficient_davies",
]

#: The Debye-Huckel constant `A` for water at 25 degC, in (mol/L)^-1/2 --
#: a commonly tabulated literature value (Atkins & de Paula, *Physical
#: Chemistry*, 11th ed., Table 5.5), derived from water's permittivity
#: and density at that temperature. `A` is temperature- and
#: solvent-dependent; this default applies specifically to aqueous
#: solutions at 25 degC.
DEBYE_HUCKEL_A_25C = 0.509


def ionic_strength(concentrations, charges) -> float:
    r"""Ionic strength :math:`I = \frac{1}{2}\sum_i c_i z_i^2`.

    Parameters
    ----------
    concentrations : array-like of float
        Molar concentration of each ion, in mol/L.
    charges : array-like of float
        Charge number of each ion (signed), same order as `concentrations`.

    Returns
    -------
    float

    Examples
    --------
    A 0.10 M solution of a 1:1 electrolyte (e.g. NaCl) has ionic strength
    equal to its concentration:

    >>> round(ionic_strength([0.10, 0.10], [1, -1]), 6)
    0.1

    A 0.10 M solution of a 1:2 electrolyte (e.g. CaCl2) has three times
    the concentration's worth of ionic strength:

    >>> round(ionic_strength([0.10, 0.20], [2, -1]), 6)
    0.3
    """
    c = np.asarray(concentrations, dtype=np.float64)
    z = np.asarray(charges, dtype=np.float64)
    return float(0.5 * np.sum(c * z**2))


def activity_coefficient_debye_huckel_limiting(z: float, I: float, A: float = DEBYE_HUCKEL_A_25C) -> float:
    r"""Debye-Huckel limiting law: :math:`\log_{10}\gamma = -A z^2 \sqrt{I}`.

    Valid only at low ionic strength (roughly :math:`I \lesssim 0.01`
    mol/L), where ion-ion interactions are dominated by long-range
    electrostatics (Debye & Huckel, 1923; Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 5.9, eq. 5.69).

    Parameters
    ----------
    z : float
        Charge number of the ion (its sign does not matter -- only
        :math:`z^2` enters).
    I : float
        Ionic strength, in mol/L.
    A : float, default :data:`DEBYE_HUCKEL_A_25C`
        Debye-Huckel constant for the solvent/temperature.

    Returns
    -------
    float
        Activity coefficient :math:`\gamma` (dimensionless, :math:`\leq 1`).

    Examples
    --------
    At zero ionic strength the activity coefficient is exactly 1 (ideal,
    infinitely dilute behavior):

    >>> round(activity_coefficient_debye_huckel_limiting(z=1, I=0.0), 6)
    1.0

    A more highly charged ion deviates further from ideality at the same
    ionic strength:

    >>> gamma_1 = activity_coefficient_debye_huckel_limiting(z=1, I=0.01)
    >>> gamma_2 = activity_coefficient_debye_huckel_limiting(z=2, I=0.01)
    >>> gamma_2 < gamma_1 < 1.0
    True
    """
    log_gamma = -A * z**2 * np.sqrt(I)
    return float(10.0**log_gamma)


def activity_coefficient_debye_huckel_extended(z: float, I: float, A: float = DEBYE_HUCKEL_A_25C, Ba: float = 1.0) -> float:
    r"""Extended Debye-Huckel law: :math:`\log_{10}\gamma = -A z^2 \sqrt{I}/(1+Ba\sqrt{I})`.

    Extends the limiting law's validity to moderate ionic strength
    (roughly :math:`I \lesssim 0.1` mol/L) by accounting for the ion's
    finite size via the dimensionless product `Ba` (the Debye-Huckel
    parameter `B` times an effective ion-size parameter `a`); `Ba=1` is a
    commonly used simplification when a specific ion-size parameter isn't
    available (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    5.9, eq. 5.70) -- this is itself still an approximation, one step
    better than the limiting law, not an exact result.

    Parameters
    ----------
    z : float
        Charge number of the ion.
    I : float
        Ionic strength, in mol/L.
    A : float, default :data:`DEBYE_HUCKEL_A_25C`
        Debye-Huckel constant for the solvent/temperature.
    Ba : float, default 1.0
        Dimensionless ion-size product; see above.

    Returns
    -------
    float
        Activity coefficient :math:`\gamma`.

    Examples
    --------
    Reduces to the limiting law as :math:`I \to 0`:

    >>> import numpy as np
    >>> I = 1e-8
    >>> extended = activity_coefficient_debye_huckel_extended(z=1, I=I)
    >>> limiting = activity_coefficient_debye_huckel_limiting(z=1, I=I)
    >>> bool(np.isclose(extended, limiting, rtol=1e-3))
    True

    At higher ionic strength, the extended law predicts a less severe
    deviation from ideality than the (over-extrapolated) limiting law:

    >>> extended_at_01 = activity_coefficient_debye_huckel_extended(z=1, I=0.1)
    >>> limiting_at_01 = activity_coefficient_debye_huckel_limiting(z=1, I=0.1)
    >>> extended_at_01 > limiting_at_01
    True
    """
    log_gamma = -A * z**2 * np.sqrt(I) / (1.0 + Ba * np.sqrt(I))
    return float(10.0**log_gamma)


def activity_coefficient_davies(z: float, I: float, A: float = DEBYE_HUCKEL_A_25C, b: float = 0.3) -> float:
    r"""Davies equation: :math:`\log_{10}\gamma = -A z^2\left(\frac{\sqrt{I}}{1+\sqrt{I}} - bI\right)`.

    An empirical extension of the Guntelberg form
    (:func:`activity_coefficient_debye_huckel_extended` with ``Ba = 1``)
    by a term linear in ionic strength, which makes :math:`\log\gamma`
    pass through a minimum and turn back up at high ionic strength, as
    measured activity coefficients do; useful up to roughly
    :math:`I \approx 0.5` mol/L with no ion-specific parameter at all.
    Davies originally used :math:`b = 0.2` (C. W. Davies, *J. Chem. Soc.*
    1938, 2093); the now-standard :math:`b = 0.3` is from his later
    monograph (C. W. Davies, *Ion Association*, Butterworths, 1962).

    Parameters
    ----------
    z : float
        Charge number of the ion.
    I : float
        Ionic strength, in mol/L.
    A : float, default :data:`DEBYE_HUCKEL_A_25C`
        Debye-Huckel constant for the solvent/temperature.
    b : float, default 0.3
        Coefficient of the linear correction term.

    Returns
    -------
    float
        Activity coefficient :math:`\gamma`.

    Examples
    --------
    With ``b = 0`` it is exactly the Guntelberg (``Ba = 1``) extended law:

    >>> g_davies = activity_coefficient_davies(z=1, I=0.05, b=0.0)
    >>> g_ext = activity_coefficient_debye_huckel_extended(z=1, I=0.05, Ba=1.0)
    >>> round(g_davies - g_ext, 12)
    0.0

    The linear term raises :math:`\gamma` again at high ionic strength:

    >>> activity_coefficient_davies(z=1, I=1.0) > activity_coefficient_davies(z=1, I=0.4)
    True
    """
    sqrt_I = np.sqrt(I)
    log_gamma = -A * z**2 * (sqrt_I / (1.0 + sqrt_I) - b * I)
    return float(10.0**log_gamma)
