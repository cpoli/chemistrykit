r"""
Laue's interference: sharp diffraction spots from a periodic array
====================================================================

Von Laue (1912) predicted that a crystal, as a periodic array of
scatterers, should diffract X-rays like a grating, and Friedrich and
Knipping saw the discrete spots. Summing the waves from :math:`N`
periodically spaced scatterers with
:func:`~chemistrykit.crystal.systems.xrd.structure_factor` gives the Laue
interference function: intensity :math:`N^2` at integer :math:`h`
(the Laue condition) and almost nothing in between, ever sharper as
:math:`N` grows.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.xrd import structure_factor

h = np.linspace(-0.5, 2.5, 1201)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for N in (3, 6, 12):
    row = [(n, 0.0, 0.0) for n in range(N)]
    intensity = np.array([abs(structure_factor((hh, 0, 0), row)) ** 2 for hh in h])
    axes[0].plot(h, intensity / N**2, label=f"N = {N}")
    peak = abs(structure_factor((1, 0, 0), row)) ** 2
    between = abs(structure_factor((0.5, 0, 0), row)) ** 2
    print(f"N={N:3d}: I(h=1) = {peak:7.1f} (= N^2 = {N * N}),  I(h=1/2) = {between:.3f}")
    assert np.isclose(peak, N * N)
axes[0].set_xlabel("h (scattering vector, reciprocal-lattice units)")
axes[0].set_ylabel(r"$|F|^2 / N^2$")
axes[0].set_title("1D Laue interference function")
axes[0].legend()

# %%
# In two dimensions the same interference gives a grid of sharp spots --
# the kind of pattern Friedrich and Knipping recorded on their plate:
N = 8
grid = [(i, j, 0.0) for i in range(N) for j in range(N)]
hk = np.linspace(-2.2, 2.2, 221)
pattern = np.array([[abs(structure_factor((hh, kk, 0), grid)) ** 2 for hh in hk] for kk in hk])
axes[1].imshow(np.log10(pattern + 1.0), extent=(hk[0], hk[-1], hk[0], hk[-1]), origin="lower", cmap="magma")
axes[1].set_xlabel("h")
axes[1].set_ylabel("k")
axes[1].set_title(f"2D array of {N}x{N} scatterers: Laue spots")
plt.tight_layout()
plt.show()
