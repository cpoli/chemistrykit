r"""
Doebereiner's platinum and Berzelius's catalysis: a catalyst that is not consumed
===================================================================================

In 1823 Doebereiner found that spongy platinum ignites a jet of hydrogen
in air at room temperature and comes out of the reaction unchanged; in
1835 Berzelius named the phenomenon *catalysis*. This example makes both
halves of that observation quantitative with illustrative (not measured)
parameters for :math:`2H_2 + O_2 \to 2H_2O`. Lowering the activation
energy on platinum
(:func:`~chemistrykit.surface.systems.catalysis.compare_catalyzed_rate`)
turns a reaction that is frozen at room temperature into a fast one, and
the :func:`~chemistrykit.surface.systems.catalysis.turnover_number` keeps
growing while the amount of platinum stays exactly the same. The
catalyst is used again and again, never used up.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.systems.catalysis import compare_catalyzed_rate, turnover_number

T = 293.15  # room temperature, K
comparison = compare_catalyzed_rate(Ea_uncatalyzed=170e3, Ea_catalyzed=60e3, T=T, A_uncatalyzed=1e11, A_catalyzed=1e9)
print(f"Uncatalyzed k at {T} K: {comparison.k_uncatalyzed:.3e} 1/s")
print(f"On platinum:            {comparison.k_catalyzed:.3e} 1/s")
print(f"Rate enhancement:       {comparison.rate_enhancement:.3e}")

# %%
# A batch of hydrogen (first order in H2 for simplicity) over a fixed
# amount of platinum. The platinum is *not* a reactant: its amount stays
# constant, while the number of H2 molecules each Pt site has converted
# (the turnover number) keeps rising.
n_H2_0 = 1.0e-2  # mol
n_Pt = 1.0e-6  # mol of surface Pt sites
t = np.linspace(0.0, 5.0 / comparison.k_catalyzed, 300)
n_H2_cat = n_H2_0 * np.exp(-comparison.k_catalyzed * t)
n_H2_uncat = n_H2_0 * np.exp(-comparison.k_uncatalyzed * t)
ton = np.array([turnover_number(n_H2_0 - n, n_Pt) for n in n_H2_cat])
print(f"\nTurnover number when the H2 is used up: {ton[-1]:.0f} cycles per Pt site")
print(f"Platinum left afterwards: {n_Pt:.1e} mol (unchanged)")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(t, n_H2_uncat / n_H2_0, label="no catalyst")
axes[0].plot(t, n_H2_cat / n_H2_0, label="spongy platinum")
axes[0].set_xlabel("time (s)")
axes[0].set_ylabel("H$_2$ remaining (fraction)")
axes[0].set_title("Doebereiner: H$_2$ + O$_2$ at room temperature")
axes[0].legend()

axes[1].plot(t, ton, color="C1", label="turnover number (cycles per Pt site)")
axes[1].set_xlabel("time (s)")
axes[1].set_ylabel("turnover number")
ax2 = axes[1].twinx()
ax2.plot(t, np.full_like(t, n_Pt * 1e6), "k--", label="Pt present (umol)")
ax2.set_ylabel("Pt present (umol)")
ax2.set_ylim(0, 2)
axes[1].set_title("Berzelius: the catalyst is not consumed")
axes[1].legend(loc="center right")
plt.tight_layout()
plt.show()
