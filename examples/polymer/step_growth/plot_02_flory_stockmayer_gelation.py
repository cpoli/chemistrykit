r"""
Flory-Stockmayer gelation: the gel point of a branching polymerization
=========================================================================

With trifunctional or higher monomers, step-growth polymerization builds
branched molecules that link into a sample-spanning network at a sharp
gel point. For self-condensation of an :math:`A_f` monomer, Flory (1941)
and Stockmayer (1943) found the weight-average degree of polymerization
:math:`\bar X_w=(1+p)/(1-(f-1)p)`
(:func:`~chemistrykit.polymer.systems.gelation.branching_weight_average_DP`)
diverges at :math:`p_c=1/(f-1)`
(:func:`~chemistrykit.polymer.systems.gelation.flory_stockmayer_gel_point`)
while the number average
(:func:`~chemistrykit.polymer.systems.gelation.branching_number_average_DP`)
stays finite. Carothers's criterion :math:`\bar X_n\to\infty`,
:math:`p_c=2/f` (:func:`~chemistrykit.polymer.systems.gelation.carothers_gel_point`),
overestimates the gel point; experiments fall between the two.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.gelation import (
    branching_number_average_DP,
    branching_weight_average_DP,
    carothers_gel_point,
    flory_stockmayer_gel_point,
)

for f in (2, 3, 4, 6):
    pc = flory_stockmayer_gel_point(f)
    Xn_gel = branching_number_average_DP(min(pc, 0.999), f)
    print(f"f = {f}: Flory-Stockmayer p_c = {pc:.3f}, Carothers p_c = {carothers_gel_point(f):.3f}, Xn at p_c = {Xn_gel:.2f}")

# %%
f = 3
pc = flory_stockmayer_gel_point(f)
p = np.linspace(0, pc - 1e-3, 400)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].semilogy(p, branching_weight_average_DP(p, f), label=r"$\bar X_w$ (diverges)")
axes[0].semilogy(p, branching_number_average_DP(p, f), label=r"$\bar X_n$ (finite)")
axes[0].axvline(pc, color="k", ls="--", lw=0.8, label=f"Flory-Stockmayer $p_c$ = {pc:.2f}")
axes[0].axvline(carothers_gel_point(f), color="crimson", ls=":", label=f"Carothers $p_c$ = {carothers_gel_point(f):.2f}")
axes[0].set_xlabel("extent of reaction p")
axes[0].set_ylabel("degree of polymerization")
axes[0].set_title("A$_3$ self-condensation")
axes[0].legend()

fs = np.linspace(2.05, 8, 100)
axes[1].plot(fs, flory_stockmayer_gel_point(fs), label="Flory-Stockmayer 1/(f-1)")
axes[1].plot(fs, carothers_gel_point(fs), "--", label="Carothers 2/f")
axes[1].set_xlabel("functionality f")
axes[1].set_ylabel("gel point $p_c$")
axes[1].set_title("Gel point vs. functionality")
axes[1].legend()
plt.tight_layout()
plt.show()
