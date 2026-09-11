r"""Solubility equilibria: Ksp, molar solubility, and the common-ion effect.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6.8, or Harris,
*Quantitative Chemical Analysis*, 9th ed., Ch. 8, throughout.

For a sparingly soluble salt :math:`M_pX_q` dissolving as
:math:`M_pX_q(s) \rightleftharpoons pM^{n+} + qX^{m-}`, at molar
solubility :math:`s` (mol/L of the salt that dissolves), :math:`[M^{n+}]
= ps` and :math:`[X^{m-}] = qs`, so

.. math::

    K_{sp} = [M^{n+}]^p[X^{m-}]^q = p^p q^q\, s^{p+q}
"""

from __future__ import annotations

from chemistrykit.solutions.utils.rootfinding import find_positive_root

__all__ = [
    "ksp_from_molar_solubility",
    "molar_solubility_from_ksp",
    "molar_solubility_with_common_ion",
]


def ksp_from_molar_solubility(s: float, cation_coeff: int, anion_coeff: int) -> float:
    r"""Compute :math:`K_{sp} = p^pq^q s^{p+q}` from a measured molar solubility.

    Parameters
    ----------
    s : float
        Molar solubility of the salt, in mol/L.
    cation_coeff : int
        Stoichiometric coefficient of the cation, `p` (e.g. 1 for AgCl, 1
        for CaF2's Ca2+, 2 for Ag2CrO4's Ag+).
    anion_coeff : int
        Stoichiometric coefficient of the anion, `q` (e.g. 1 for AgCl, 2
        for CaF2's F-, 1 for Ag2CrO4's CrO4 2-).

    Returns
    -------
    float

    Examples
    --------
    AgCl (1:1 salt) with molar solubility 1.3e-5 mol/L:

    >>> round(ksp_from_molar_solubility(s=1.3e-5, cation_coeff=1, anion_coeff=1), 12)
    1.69e-10

    CaF2 (1:2 salt), Ksp = [Ca2+][F-]^2 = s*(2s)^2 = 4s^3:

    >>> round(ksp_from_molar_solubility(s=2.1e-4, cation_coeff=1, anion_coeff=2), 12)
    3.7e-11
    """
    p, q = cation_coeff, anion_coeff
    return (p**p) * (q**q) * s ** (p + q)


def molar_solubility_from_ksp(Ksp: float, cation_coeff: int, anion_coeff: int) -> float:
    r"""Invert :func:`ksp_from_molar_solubility`: :math:`s = (K_{sp}/(p^pq^q))^{1/(p+q)}`.

    Parameters
    ----------
    Ksp : float
        Solubility product.
    cation_coeff : int
        Stoichiometric coefficient of the cation, `p`.
    anion_coeff : int
        Stoichiometric coefficient of the anion, `q`.

    Returns
    -------
    float
        Molar solubility, in mol/L.

    Examples
    --------
    Round-trips with :func:`ksp_from_molar_solubility`:

    >>> Ksp = ksp_from_molar_solubility(s=1.3e-5, cation_coeff=1, anion_coeff=1)
    >>> round(molar_solubility_from_ksp(Ksp, cation_coeff=1, anion_coeff=1), 12)
    1.3e-05
    """
    p, q = cation_coeff, anion_coeff
    return (Ksp / ((p**p) * (q**q))) ** (1.0 / (p + q))


def molar_solubility_with_common_ion(Ksp: float, cation_coeff: int, anion_coeff: int, common_ion_conc: float, common_ion: str = "cation") -> float:
    r"""Molar solubility of :math:`M_pX_q` in the presence of an added common ion.

    Le Chatelier's principle predicts that adding a common ion
    (independently, e.g. from a soluble salt sharing that ion) suppresses
    the target salt's solubility relative to :func:`molar_solubility_from_ksp`
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6.8b). This
    solves the exact polynomial

    .. math::

        K_{sp} = (ps + C_0)^p (qs)^q \quad\text{(common cation, added at concentration } C_0\text{)}

    or the mirror-image equation for a common anion, via
    :func:`chemistrykit.solutions.utils.rootfinding.find_positive_root`
    (exact for any p, q -- not just the textbook 1:1-salt quadratic
    special case).

    Parameters
    ----------
    Ksp : float
        Solubility product of the salt whose solubility is being found.
    cation_coeff : int
        Stoichiometric coefficient of the cation, `p`.
    anion_coeff : int
        Stoichiometric coefficient of the anion, `q`.
    common_ion_conc : float
        Concentration of the common ion contributed by an independent
        (fully dissociated) source, in mol/L.
    common_ion : {"cation", "anion"}
        Which ion is shared with the independent source.

    Returns
    -------
    float
        Molar solubility of the salt, in mol/L.

    Examples
    --------
    AgCl (Ksp = 1.8e-10) is much less soluble in 0.10 M NaCl (common Cl-
    ion) than in pure water:

    >>> s_pure = molar_solubility_from_ksp(1.8e-10, cation_coeff=1, anion_coeff=1)
    >>> s_common = molar_solubility_with_common_ion(1.8e-10, cation_coeff=1, anion_coeff=1, common_ion_conc=0.10, common_ion="anion")
    >>> s_common < s_pure
    True
    >>> round(s_common, 12)  # dominated by the 0.10 M Cl- already present: s ~= Ksp/0.10
    1.8e-09
    """
    p, q = cation_coeff, anion_coeff

    def residual(s):
        if common_ion == "cation":
            cation = p * s + common_ion_conc
            anion = q * s
        elif common_ion == "anion":
            cation = p * s
            anion = q * s + common_ion_conc
        else:
            raise ValueError("common_ion must be 'cation' or 'anion'")
        return cation**p * anion**q - Ksp

    # An upper bound for s: the common-ion-free solubility (adding a
    # common ion can only suppress solubility further).
    s_free = molar_solubility_from_ksp(Ksp, p, q)
    return find_positive_root(residual, lo=1e-20, hi=max(s_free, 1e-10))
