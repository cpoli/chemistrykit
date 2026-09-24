r"""
The Eley-Rideal mechanism: a gas molecule strikes an adsorbed one
===================================================================

In the Eley-Rideal mechanism (Eley and Rideal, 1940) only reactant A
adsorbs. B reacts with adsorbed A straight from the gas phase, so the
rate is :math:`k\theta_AP_B`
(:func:`~chemistrykit.surface.systems.eley_rideal.er_rate`). A and B never
compete for sites, so raising :math:`P_A` can only help: the rate rises
and levels off, and it stays first order in :math:`P_B`. In the dual-site
Langmuir-Hinshelwood mechanism
(:func:`~chemistrykit.surface.systems.langmuir_hinshelwood.lh_rate_dual_site`),
too much A pushes B off the surface and the rate turns over. That
difference is the standard kinetic test for telling the two mechanisms
apart.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.eley_rideal import er_rate
from chemistrykit.surface.systems.langmuir_hinshelwood import lh_rate_dual_site

k, K_A, K_B, P_B = 4.0, 2.0, 5.0, 0.2
P_A = np.linspace(0.001, 20.0, 400)
rate_er = er_rate(k, K_A, P_A, P_B)
rate_lh = lh_rate_dual_site(k, K_A, P_A, K_B, P_B)
print(f"Eley-Rideal rate is monotonic in P_A: {bool(np.all(np.diff(rate_er) > 0))}")
print(f"Eley-Rideal saturates at k*P_B = {k * P_B:.2f}; at P_A=20: {rate_er[-1]:.3f}")
print(f"Langmuir-Hinshelwood peaks at P_A = {P_A[np.argmax(rate_lh)]:.3f}, then falls")

# %%
P_B_scan = np.linspace(0.0, 5.0, 200)
er_vs_B = er_rate(k, K_A, 1.0, P_B_scan)
lh_vs_B = lh_rate_dual_site(k, K_A, 1.0, K_B, P_B_scan)

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(P_A, rate_er / rate_er.max(), label="Eley-Rideal")
axes[0].plot(P_A, rate_lh / rate_lh.max(), label="Langmuir-Hinshelwood")
axes[0].set_xlabel("$P_A$")
axes[0].set_ylabel("rate (normalized)")
axes[0].set_title("Rate vs. pressure of adsorbed reactant A")
axes[0].legend()
axes[1].plot(P_B_scan, er_vs_B / er_vs_B.max(), label="Eley-Rideal: first order in $P_B$")
axes[1].plot(P_B_scan, lh_vs_B / lh_vs_B.max(), label="Langmuir-Hinshelwood")
axes[1].set_xlabel("$P_B$")
axes[1].set_ylabel("rate (normalized)")
axes[1].set_title("Rate vs. pressure of gas-phase reactant B")
axes[1].legend()
plt.tight_layout()
plt.show()
