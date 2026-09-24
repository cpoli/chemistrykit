r"""
Frenkel defects: vacancy-interstitial pairs in silver chloride
================================================================

Frenkel (1926) argued that at any :math:`T>0` some ions leave their sites
for interstitial gaps, because the configurational entropy gained outweighs
the enthalpy cost at low concentration. The equilibrium number is
:math:`n_F=\sqrt{NN_i}\exp(-\Delta H_F/2k_BT)`
(:func:`~chemistrykit.crystal.systems.defects.frenkel_defect_concentration`).
AgCl, where small Ag+ ions readily occupy interstitial sites, is the
textbook Frenkel-dominated solid (:math:`\Delta H_F\approx1.4` eV).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, K_B
from chemistrykit.crystal.systems.defects import frenkel_defect_concentration

N = 2.3e22  # Ag+ sites per cm^3 in AgCl (approximate)
delta_h = 1.4 * ELEMENTARY_CHARGE  # J
T = np.linspace(300.0, 700.0, 200)  # AgCl melts at 728 K

for T_check in (300.0, 500.0, 700.0):
    nf = frenkel_defect_concentration(N, 2.0 * N, delta_h, T_check)
    print(f"T = {T_check:4.0f} K: fraction of Ag+ displaced to interstitials = {nf / N:.2e}")

# %%
# The Arrhenius slope of ln(n_F) vs 1/T is -Delta H_F / (2 k_B): the factor
# 2 because the vacancy and the interstitial are created together.
n = frenkel_defect_concentration(N, 2.0 * N, delta_h, T)
slope = np.polyfit(1.0 / T, np.log(n), 1)[0]
print(f"\nfitted slope {slope:.1f} K,  -Delta H_F / 2k_B = {-delta_h / (2.0 * K_B):.1f} K")
assert np.isclose(slope, -delta_h / (2.0 * K_B))

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ratio in (0.5, 2.0, 8.0):
    axes[0].semilogy(T, frenkel_defect_concentration(N, ratio * N, delta_h, T) / N, label=f"N_i = {ratio:g} N")
axes[0].set_xlabel("T (K)")
axes[0].set_ylabel(r"$n_F/N$")
axes[0].set_title(r"Frenkel defects: more interstitial sites, more defects ($\propto\sqrt{N_i}$)")
axes[0].legend()
axes[1].plot(1000.0 / T, np.log(n / N))
axes[1].set_xlabel("1000 / T (1/K)")
axes[1].set_ylabel(r"$\ln(n_F/N)$")
axes[1].set_title(r"Arrhenius plot: slope $-\Delta H_F/2k_B$")
plt.tight_layout()
plt.show()
