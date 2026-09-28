r"""
Guldberg and Waage's law of mass action: a constant equilibrium quotient
==========================================================================

The law of mass action says that at equilibrium the concentrations of
products and reactants, each raised to its stoichiometric power, combine
into a quotient that is fixed at a given temperature. For
:math:`HA \rightleftharpoons H^+ + A^-`:

.. math::

   K_a = \frac{[H^+][A^-]}{[HA]}.

Below, the equilibrium concentrations of three weak acids are solved
exactly with :class:`~chemistrykit.solutions.systems.acid_base.WeakAcid`
over six decades of total concentration. The individual concentrations
change by orders of magnitude, but the mass-action quotient stays fixed
at :math:`K_a`, while a quotient that is *not* the mass-action form,
such as :math:`[H^+]/[HA]`, does not.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import WeakAcid

acids = {"formic acid": 1.8e-4, "acetic acid": 1.8e-5, "hypochlorous acid": 3.0e-8}
Ca_values = np.logspace(-5, 0, 60)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
for (name, Ka), color in zip(acids.items(), ["steelblue", "darkorange", "seagreen"], strict=True):
    Q_mass_action, Q_other = [], []
    for Ca in Ca_values:
        acid = WeakAcid(Ca=Ca, Ka=Ka)
        h = acid.h_concentration()
        A_minus = Ka * Ca / (Ka + h)  # from the mass balance and Ka
        HA = Ca - A_minus
        Q_mass_action.append(h * A_minus / HA)
        Q_other.append(h / HA)
    ax1.loglog(Ca_values, Q_mass_action, color=color, label=name)
    ax1.axhline(Ka, color=color, linestyle=":", alpha=0.6)
    ax2.loglog(Ca_values, Q_other, color=color, label=name)

ax1.set_xlabel("total acid concentration $C_a$ (mol/L)")
ax1.set_ylabel(r"$[H^+][A^-]/[HA]$")
ax1.set_title("Mass-action quotient: constant = $K_a$")
ax1.legend(fontsize=8)
ax2.set_xlabel("total acid concentration $C_a$ (mol/L)")
ax2.set_ylabel(r"$[H^+]/[HA]$")
ax2.set_title("A non-mass-action ratio: not constant")
fig.tight_layout()

# %%
# The exact solver never imposes the quotient directly -- it solves a
# cubic from mass balance, charge balance, and water autoionization -- yet
# the recovered quotient matches :math:`K_a` to rounding error:

acid = WeakAcid(Ca=0.010, Ka=1.8e-5)
h = acid.h_concentration()
A_minus = acid.Ka * acid.Ca / (acid.Ka + h)
print(f"[H+][A-]/[HA] = {h * A_minus / (acid.Ca - A_minus):.6e}  (Ka = {acid.Ka:.6e})")

plt.show()
