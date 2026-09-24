r"""
The Sackur-Tetrode entropy of a monatomic gas
================================================

Sackur (1911) and Tetrode (1912) obtained the absolute entropy of an ideal
monatomic gas by dividing classical phase space into cells of size
:math:`h^3`:

.. math::

    S_m = R\left[\ln\left(\frac{k_BT}{P\Lambda^3}\right)+\frac52\right],
    \qquad \Lambda=\frac{h}{\sqrt{2\pi mk_BT}}

:func:`~chemistrykit.statmech.sackur_tetrode_entropy` evaluates it from
the atomic mass, temperature, and pressure alone. Below it is compared
with the standard molar entropies of the noble gases (Atkins & de Paula,
Table 13.1 lists argon at about 154.8 J/(mol K)), and its
:math:`\frac32R\ln m` mass dependence is shown directly.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.constants import R
from chemistrykit.statmech import sackur_tetrode_entropy

# Standard molar entropies at 298.15 K and 1 bar, J/(mol K) (NIST-JANAF / CODATA key values)
noble_gases = {"He": (4.0026, 126.15), "Ne": (20.180, 146.33), "Ar": (39.948, 154.85), "Kr": (83.798, 164.09), "Xe": (131.29, 169.69)}

names = list(noble_gases)
masses = np.array([noble_gases[n][0] for n in names])
S_exp = np.array([noble_gases[n][1] for n in names])
S_st = np.array([sackur_tetrode_entropy(mass=m * sc.atomic_mass, T=298.15, P=1.0e5) for m in masses])

for name, s_calc, s_exp in zip(names, S_st, S_exp):
    print(f"{name:2s}: Sackur-Tetrode {s_calc:7.2f}  tabulated {s_exp:7.2f} J/(mol K)")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
m_grid = np.logspace(0.3, 2.3, 200)
axes[0].semilogx(m_grid, [sackur_tetrode_entropy(mass=m * sc.atomic_mass, T=298.15, P=1.0e5) for m in m_grid], color="steelblue", label="Sackur-Tetrode")
axes[0].semilogx(masses, S_exp, "o", color="crimson", label="tabulated (298.15 K, 1 bar)")
for name, m, s in zip(names, masses, S_exp):
    axes[0].annotate(name, (m, s), textcoords="offset points", xytext=(5, -12))
axes[0].set_xlabel("atomic mass (u)")
axes[0].set_ylabel(r"$S_m$ (J mol$^{-1}$ K$^{-1}$)")
axes[0].set_title(r"Noble-gas entropies: slope $\frac{3}{2}R$ per $\ln m$")
axes[0].legend()

T = np.linspace(50.0, 1500.0, 200)
argon = 39.948 * sc.atomic_mass
S_T = np.array([sackur_tetrode_entropy(mass=argon, T=t, P=1.0e5) for t in T])
axes[1].plot(T, S_T, color="steelblue", label="argon, 1 bar")
axes[1].plot(T, S_T[0] + 2.5 * R * np.log(T / T[0]), color="gray", linestyle="--", label=r"$\frac{5}{2}R\ln T$ at constant $P$")
axes[1].set_xlabel("T (K)")
axes[1].set_ylabel(r"$S_m$ (J mol$^{-1}$ K$^{-1}$)")
axes[1].set_title("Temperature dependence")
axes[1].legend()
fig.tight_layout()

plt.show()
