r"""
Guntelberg's extended law: correcting Debye-Huckel for finite ion size
========================================================================

The limiting law treats ions as point charges and over-corrects above
:math:`I \approx 0.01` mol/L. Guntelberg's extended form adds a
finite-size denominator with a single universal size parameter,
:math:`Ba \approx 1`,

.. math::

   \log_{10}\gamma = \frac{-A z^2 \sqrt{I}}{1 + \sqrt{I}},

implemented by
:func:`~chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_extended`.
It coincides with the limiting law
(:func:`~chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_limiting`)
at high dilution and stays much closer to reality an order of magnitude
further in ionic strength.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.activity import (
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
)

I_range = np.logspace(-4, -0.3, 100)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
for z, color in [(1, "steelblue"), (2, "darkorange"), (3, "seagreen")]:
    g_lim = np.array([activity_coefficient_debye_huckel_limiting(z, I) for I in I_range])
    g_gun = np.array([activity_coefficient_debye_huckel_extended(z, I, Ba=1.0) for I in I_range])
    ax1.semilogx(I_range, g_lim, linestyle="--", color=color, alpha=0.6, label=f"z={z}, limiting")
    ax1.semilogx(I_range, g_gun, color=color, label=f"z={z}, Guntelberg")
    ax2.semilogx(I_range, g_gun / g_lim, color=color, label=f"z={z}")
ax1.set_xlabel("ionic strength I (mol/L)")
ax1.set_ylabel(r"activity coefficient $\gamma$")
ax1.set_title("Limiting law vs. Guntelberg's extended law")
ax1.legend(fontsize=8)
ax2.axvline(0.01, color="gray", linestyle=":", label="I = 0.01 M")
ax2.set_xlabel("ionic strength I (mol/L)")
ax2.set_ylabel(r"$\gamma_{\mathrm{Guntelberg}}/\gamma_{\mathrm{limiting}}$")
ax2.set_title("Size of the finite-ion-size correction")
ax2.legend(fontsize=8)
fig.tight_layout()

# %%
# Effect of the ion-size parameter Ba for a singly charged ion at I = 0.1 M
# (Guntelberg's choice is Ba = 1):

for Ba in (0.0, 0.5, 1.0, 1.5, 2.0):
    print(f"Ba = {Ba:.1f}: gamma = {activity_coefficient_debye_huckel_extended(1, 0.1, Ba=Ba):.4f}")

plt.show()
