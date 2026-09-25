r"""
Ostwald's dilution law: weak acids dissociate more when diluted
=================================================================

Dissolving less acid in more water increases the *fraction* of it that
dissociates, even though the absolute :math:`[H^+]` goes down -- the
dilution law Wilhelm Ostwald derived in 1888 by applying the mass-action
law to Arrhenius's dissociation equilibrium,

.. math::

   K_a = \frac{\alpha^2 C_a}{1-\alpha}.

Here the degree of dissociation :math:`\alpha` of acetic acid is
computed from the exact cubic charge-balance solution in
:class:`~chemistrykit.solutions.systems.acid_base.WeakAcid` and compared
against Ostwald's own closed-form solution of the quadratic above (which
ignores water autoionization).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import WeakAcid

Ka = 1.8e-5  # acetic acid
Ca_values = np.logspace(-6, 0, 80)
alpha_exact = np.array([WeakAcid(Ca=Ca, Ka=Ka).percent_dissociation() / 100.0 for Ca in Ca_values])

# Ostwald's quadratic alpha^2 Ca + Ka alpha - Ka = 0, positive root:
alpha_ostwald = (-Ka + np.sqrt(Ka**2 + 4.0 * Ka * Ca_values)) / (2.0 * Ca_values)

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(Ca_values, 100.0 * alpha_exact, color="steelblue", label="exact (WeakAcid)")
ax.semilogx(Ca_values, 100.0 * alpha_ostwald, color="darkorange", linestyle="--", label=r"Ostwald: $K_a = \alpha^2 C_a/(1-\alpha)$")
ax.set_xlabel("total acid concentration $C_a$ (mol/L)")
ax.set_ylabel(r"percent dissociated, $100\alpha$")
ax.set_title("Ostwald dilution law: acetic acid")
ax.legend()
fig.tight_layout()

# %%
# The "Ostwald constant" :math:`\alpha^2 C_a/(1-\alpha)` computed from the
# exact degree of dissociation stays equal to :math:`K_a` across four
# decades of dilution -- the experimental test Ostwald made with
# conductivity data for dozens of weak acids:

for Ca, alpha in zip(Ca_values[::20], alpha_exact[::20], strict=True):
    print(f"Ca = {Ca:.1e} M: alpha = {alpha:.4f}, alpha^2 Ca/(1-alpha) = {alpha**2 * Ca / (1 - alpha):.3e}")

plt.show()
