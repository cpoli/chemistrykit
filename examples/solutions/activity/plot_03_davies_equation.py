r"""
The Davies equation: activity coefficients that turn back up
==============================================================

Measured activity coefficients fall with ionic strength, pass through a
minimum, and rise again in concentrated solutions -- behaviour no pure
screening law reproduces. Davies added an empirical term linear in
:math:`I` to Guntelberg's form,

.. math::

   \log_{10}\gamma = -A z^2\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3\,I\right),

implemented by
:func:`~chemistrykit.solutions.systems.activity.activity_coefficient_davies`.
The comparison below shows the minimum near :math:`I \approx 0.4` mol/L
and the upturn that the Debye-Huckel family misses.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.activity import (
    activity_coefficient_davies,
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
)

I_range = np.linspace(1e-4, 2.0, 400)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(I_range, [activity_coefficient_debye_huckel_limiting(1, I) for I in I_range], color="gray", linestyle=":", label="Debye-Huckel limiting")
ax.plot(I_range, [activity_coefficient_debye_huckel_extended(1, I) for I in I_range], color="darkorange", linestyle="--", label="Guntelberg (Ba = 1)")
ax.plot(I_range, [activity_coefficient_davies(1, I, b=0.2) for I in I_range], color="seagreen", alpha=0.7, label="Davies, b = 0.2 (1938)")
ax.plot(I_range, [activity_coefficient_davies(1, I) for I in I_range], color="steelblue", label="Davies, b = 0.3 (1962)")
ax.set_xlabel("ionic strength I (mol/L)")
ax.set_ylabel(r"activity coefficient $\gamma$, z = 1")
ax.set_title("Davies equation vs. the Debye-Huckel family")
ax.legend(fontsize=8)
fig.tight_layout()

# %%
# Location of the Davies minimum for b = 0.3, found on a fine grid:

I_fine = np.linspace(0.01, 2.0, 20001)
gamma_fine = np.array([activity_coefficient_davies(1, I) for I in I_fine])
i_min = int(np.argmin(gamma_fine))
print(f"Davies minimum: gamma = {gamma_fine[i_min]:.4f} at I = {I_fine[i_min]:.3f} mol/L")

plt.show()
