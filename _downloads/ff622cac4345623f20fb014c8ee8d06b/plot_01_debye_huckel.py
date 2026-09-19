r"""
Debye-Huckel limiting vs. extended activity-coefficient laws
===============================================================

Comparing the Debye-Huckel limiting law
(:func:`~chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_limiting`)
against the extended law
(:func:`~chemistrykit.solutions.systems.activity.activity_coefficient_debye_huckel_extended`)
for ions of different charge, across a range of ionic strength. The two
laws agree closely at low ionic strength but diverge above roughly 0.01
mol/L, where the limiting law over-corrects.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.activity import (
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
    ionic_strength,
)

I_range = np.logspace(-4, -0.3, 100)

fig, ax = plt.subplots(figsize=(7, 5))
for z, color in [(1, "steelblue"), (2, "darkorange"), (3, "seagreen")]:
    gamma_limiting = [activity_coefficient_debye_huckel_limiting(z, I) for I in I_range]
    gamma_extended = [activity_coefficient_debye_huckel_extended(z, I) for I in I_range]
    ax.semilogx(I_range, gamma_limiting, linestyle="--", color=color, alpha=0.6, label=f"z={z}, limiting")
    ax.semilogx(I_range, gamma_extended, color=color, label=f"z={z}, extended")

ax.set_xlabel("ionic strength I (mol/L)")
ax.set_ylabel("activity coefficient")
ax.set_title("Debye-Huckel limiting vs. extended law")
ax.legend(fontsize=8)
fig.tight_layout()

# %%
# The ionic strength of a real mixed-electrolyte solution, e.g. 0.05 M
# CaCl2 plus 0.02 M NaCl:

I_mixed = ionic_strength(concentrations=[0.05, 0.10, 0.02, 0.02], charges=[2, -1, 1, -1])
print(f"Ionic strength of 0.05 M CaCl2 + 0.02 M NaCl: {I_mixed:.4f} mol/L")

plt.show()
