r"""
Kratky-Porod worm-like chain: from rigid rod to random coil
==============================================================

A stiff polymer (DNA, many polysaccharides) bends smoothly rather than at
free joints. Kratky and Porod's (1949) worm-like chain, with persistence
length :math:`P`, gives
(:func:`~chemistrykit.polymer.systems.chain_statistics.worm_like_chain_mean_square_end_to_end`)

.. math::

    \langle R^2\rangle = 2PL\left[1-\frac{P}{L}\left(1-e^{-L/P}\right)\right]

a rigid rod (:math:`R^2=L^2`) for contour lengths :math:`L\ll P` and an
ideal Kuhn coil (:math:`R^2=2PL`, Kuhn length :math:`b=2P`) for
:math:`L\gg P`. Double-stranded DNA, with :math:`P\approx50` nm, is used
as an illustration.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.chain_statistics import worm_like_chain_mean_square_end_to_end

P = 50.0  # nm, dsDNA persistence length
L = np.logspace(0, 5, 200)  # contour length in nm
R2 = worm_like_chain_mean_square_end_to_end(L, P)

for Lv in (10.0, 50.0, 500.0, 50000.0):
    r = float(worm_like_chain_mean_square_end_to_end(Lv, P))
    print(f"L = {Lv:8.0f} nm:  R = {np.sqrt(r):8.1f} nm   R^2/L^2 = {r / Lv**2:.3f}   R^2/(2PL) = {r / (2 * P * Lv):.3f}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].loglog(L, np.sqrt(R2), label="worm-like chain")
axes[0].loglog(L, L, "k--", lw=0.8, label="rigid rod, R = L")
axes[0].loglog(L, np.sqrt(2 * P * L), "k:", lw=0.8, label=r"Kuhn coil, $R=\sqrt{2PL}$")
axes[0].axvline(P, color="gray", lw=0.6)
axes[0].set_xlabel("contour length L (nm)")
axes[0].set_ylabel(r"$\sqrt{\langle R^2\rangle}$ (nm)")
axes[0].set_title("dsDNA, P = 50 nm")
axes[0].legend()

for Pv in (5.0, 50.0, 500.0):
    axes[1].semilogx(L / Pv, worm_like_chain_mean_square_end_to_end(L, Pv) / L**2, label=f"P = {Pv:g} nm")
axes[1].set_xlabel("L / P")
axes[1].set_ylabel(r"$\langle R^2\rangle / L^2$")
axes[1].set_title("Stiffness collapses onto one curve in L/P")
axes[1].legend()
plt.tight_layout()
plt.show()
