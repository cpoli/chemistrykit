r"""Theories of the rate constant itself: collisions, diffusion, and transition states.

Three classic routes from molecular properties to a bimolecular or
unimolecular rate constant, each a closed-form formula:

* :func:`collision_theory_rate_constant` -- hard-sphere collision theory
  (M. Trautz, *Z. Anorg. Allg. Chem.* 96, 1 (1916); W. C. McC. Lewis,
  *J. Chem. Soc., Trans.* 113, 471 (1918)).
* :func:`smoluchowski_rate_constant`,
  :func:`smoluchowski_transient_rate_constant`, and
  :func:`diffusion_limited_rate_constant` -- the diffusion-controlled
  encounter rate (M. v. Smoluchowski, *Z. Phys. Chem.* 92, 129 (1917)).
* :func:`eyring_rate_constant` and :func:`fit_eyring` -- the
  thermodynamic (Eyring-Polanyi) form of transition-state theory (H.
  Eyring, *J. Chem. Phys.* 3, 107 (1935)).

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 18 ("Reaction
dynamics") for all three.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import K_B, NA, H, R
from chemistrykit.kinetics.utils.linear_regression import linear_fit

__all__ = [
    "collision_theory_rate_constant",
    "smoluchowski_rate_constant",
    "smoluchowski_transient_rate_constant",
    "diffusion_limited_rate_constant",
    "eyring_rate_constant",
    "EyringFit",
    "fit_eyring",
]


def collision_theory_rate_constant(T, sigma: float, reduced_mass: float, Ea: float = 0.0, steric_factor: float = 1.0):
    r"""Hard-sphere collision-theory rate constant of a bimolecular reaction.

    Every collision between reactant molecules with relative kinetic
    energy along the line of centres above :math:`E_a` reacts (times an
    empirical steric factor :math:`P`):

    .. math::

        k = P\,\sigma \sqrt{\frac{8 k_B T}{\pi \mu}}\; N_A\, e^{-E_a/RT}

    where :math:`\sigma = \pi d^2` is the collision cross-section and
    :math:`\sqrt{8k_BT/\pi\mu}` the mean relative speed (Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Topic 18A).

    Parameters
    ----------
    T : float or array-like of float
        Absolute temperature(s), in K.
    sigma : float
        Collision cross-section, in m^2.
    reduced_mass : float
        Reduced mass :math:`\mu = m_A m_B/(m_A + m_B)` of the colliding
        pair, in kg.
    Ea : float, default 0.0
        Activation (threshold) energy, in J/mol.
    steric_factor : float, default 1.0
        Steric factor :math:`P`.

    Returns
    -------
    float or ndarray
        Rate constant, in m^3 mol^-1 s^-1.

    Examples
    --------
    With no barrier, :math:`k \propto \sqrt{T}` exactly:

    >>> import numpy as np
    >>> k = collision_theory_rate_constant(np.array([300.0, 1200.0]), sigma=4e-19, reduced_mass=1e-26)
    >>> round(float(k[1] / k[0]), 12)
    2.0
    """
    T = np.asarray(T, dtype=np.float64)
    mean_speed = np.sqrt(8.0 * K_B * T / (np.pi * reduced_mass))
    return steric_factor * sigma * mean_speed * NA * np.exp(-Ea / (R * T))


def smoluchowski_rate_constant(D: float, R_contact: float):
    r"""Smoluchowski's steady-state diffusion-limited encounter rate constant.

    Two species diffusing with combined diffusion coefficient
    :math:`D = D_A + D_B` react on first contact at separation
    :math:`R^* = R_A + R_B`; solving the steady-state diffusion equation
    around one reactant gives

    .. math::

        k_D = 4\pi D R^* N_A

    Parameters
    ----------
    D : float
        Sum of the two reactants' diffusion coefficients, in m^2/s.
    R_contact : float
        Reaction (contact) distance, in m.

    Returns
    -------
    float
        Rate constant, in m^3 mol^-1 s^-1.

    Examples
    --------
    Two small molecules in water (D ~ 2e-9 m^2/s each, contact at 0.5 nm)
    react at about 1.5e10 L mol^-1 s^-1 (1 m^3 = 1000 L):

    >>> k = smoluchowski_rate_constant(D=4e-9, R_contact=5e-10)
    >>> round(k * 1000 / 1e10, 2)
    1.51
    """
    return 4.0 * np.pi * D * R_contact * NA


def smoluchowski_transient_rate_constant(D: float, R_contact: float, t):
    r"""Smoluchowski's time-dependent diffusion-limited rate constant.

    Starting from a uniform (random) distribution of reactants, the
    depletion zone around each reactant has not yet formed, so the rate
    coefficient starts higher than :math:`k_D` and decays towards it:

    .. math::

        k(t) = 4\pi D R^* N_A \left(1 + \frac{R^*}{\sqrt{\pi D t}}\right)

    Parameters
    ----------
    D : float
        Sum of the two reactants' diffusion coefficients, in m^2/s.
    R_contact : float
        Reaction (contact) distance, in m.
    t : float or array-like of float
        Time(s) since mixing, in s (must be positive).

    Returns
    -------
    float or ndarray
        Rate constant(s), in m^3 mol^-1 s^-1.

    Examples
    --------
    At :math:`t = R^{*2}/(\pi D)` the transient term exactly doubles the
    steady-state rate:

    >>> import numpy as np
    >>> D, Rc = 4e-9, 5e-10
    >>> ratio = smoluchowski_transient_rate_constant(D, Rc, Rc**2 / (np.pi * D)) / smoluchowski_rate_constant(D, Rc)
    >>> round(float(ratio), 12)
    2.0
    """
    t = np.asarray(t, dtype=np.float64)
    return smoluchowski_rate_constant(D, R_contact) * (1.0 + R_contact / np.sqrt(np.pi * D * t))


def diffusion_limited_rate_constant(T, viscosity: float):
    r"""Diffusion-limited rate constant from the solvent viscosity alone.

    Combining :func:`smoluchowski_rate_constant` with the Stokes-Einstein
    relation :math:`D = k_B T/(6\pi\eta r)` for two reactants of equal
    radius :math:`r` (so :math:`D = D_A + D_B = 2k_BT/(6\pi\eta r)` and
    :math:`R^* = 2r`) makes the radius cancel:

    .. math::

        k_D = \frac{8RT}{3\eta}

    Parameters
    ----------
    T : float or array-like of float
        Absolute temperature(s), in K.
    viscosity : float
        Solvent dynamic viscosity :math:`\eta`, in Pa s.

    Returns
    -------
    float or ndarray
        Rate constant, in m^3 mol^-1 s^-1.

    Examples
    --------
    Water at 298 K (:math:`\eta \approx 8.9\times 10^{-4}` Pa s) gives
    the familiar :math:`\sim 7\times 10^{9}` L mol^-1 s^-1:

    >>> k = diffusion_limited_rate_constant(298.15, 8.9e-4)
    >>> round(float(k) * 1000 / 1e9, 2)
    7.43
    """
    T = np.asarray(T, dtype=np.float64)
    return 8.0 * R * T / (3.0 * viscosity)


def eyring_rate_constant(T, dH: float, dS: float, kappa: float = 1.0):
    r"""Eyring (transition-state theory) rate constant.

    .. math::

        k = \kappa \frac{k_B T}{h}\, e^{\Delta S^{\ddagger}/R}\, e^{-\Delta H^{\ddagger}/RT}

    i.e. :math:`k = \kappa (k_BT/h) e^{-\Delta G^{\ddagger}/RT}` with
    :math:`\Delta G^{\ddagger} = \Delta H^{\ddagger} - T\Delta S^{\ddagger}`
    (H. Eyring, *J. Chem. Phys.* 3, 107 (1935); Atkins & de Paula,
    *Physical Chemistry*, 11th ed., Topic 18C). Written for a unimolecular
    step (or with the standard-state concentration factor absorbed into
    :math:`\Delta S^{\ddagger}`), so `k` is in s^-1.

    Parameters
    ----------
    T : float or array-like of float
        Absolute temperature(s), in K.
    dH : float
        Enthalpy of activation :math:`\Delta H^{\ddagger}`, in J/mol.
    dS : float
        Entropy of activation :math:`\Delta S^{\ddagger}`, in J mol^-1 K^-1.
    kappa : float, default 1.0
        Transmission coefficient.

    Returns
    -------
    float or ndarray
        Rate constant(s), in s^-1.

    Examples
    --------
    With zero activation enthalpy and entropy, `k` is the universal
    frequency :math:`k_BT/h` (about 6.25e12 s^-1 at 300 K):

    >>> round(float(eyring_rate_constant(300.0, dH=0.0, dS=0.0)) / 1e12, 3)
    6.251
    """
    T = np.asarray(T, dtype=np.float64)
    return kappa * K_B * T / H * np.exp(dS / R) * np.exp(-dH / (R * T))


@dataclass
class EyringFit:
    """Result of fitting rate-constant-vs-temperature data to the Eyring equation."""

    dH: float
    """float: Fitted enthalpy of activation, in J/mol."""

    dS: float
    """float: Fitted entropy of activation, in J mol^-1 K^-1."""

    r_squared: float
    """float: Coefficient of determination of the linear (ln(k/T) vs 1/T) fit."""

    def dG(self, T):
        """Gibbs energy of activation ``dH - T*dS`` at temperature(s) `T`, in J/mol.

        Parameters
        ----------
        T : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return self.dH - np.asarray(T, dtype=np.float64) * self.dS

    def predict(self, T):
        """Evaluate the fitted Eyring equation at temperature(s) `T`.

        Parameters
        ----------
        T : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return eyring_rate_constant(T, self.dH, self.dS)


def fit_eyring(T, k) -> EyringFit:
    r"""Fit rate constant vs. temperature data to the Eyring equation.

    Linearizes :math:`\ln(k/T) = \ln(k_B/h) + \Delta S^{\ddagger}/R -
    (\Delta H^{\ddagger}/R)(1/T)` and fits by ordinary least squares --
    the "Eyring plot": slope :math:`-\Delta H^{\ddagger}/R`, intercept
    :math:`\ln(k_B/h) + \Delta S^{\ddagger}/R`.

    Parameters
    ----------
    T : array-like of float
        Absolute temperatures, in K (at least 2 distinct values).
    k : array-like of float
        First-order rate constants at each temperature, in s^-1.

    Returns
    -------
    EyringFit

    Examples
    --------
    >>> import numpy as np
    >>> T = np.linspace(280.0, 360.0, 5)
    >>> fit = fit_eyring(T, eyring_rate_constant(T, dH=80e3, dS=-20.0))
    >>> round(fit.dH, 3), round(fit.dS, 6)
    (80000.0, -20.0)
    """
    T = np.asarray(T, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    fit = linear_fit(1.0 / T, np.log(k / T))
    dH = -fit.slope * R
    dS = (fit.intercept - np.log(K_B / H)) * R
    return EyringFit(dH=float(dH), dS=float(dS), r_squared=fit.r_squared)
