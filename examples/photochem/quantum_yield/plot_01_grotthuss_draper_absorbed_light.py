r"""
Grotthuss-Draper law: only absorbed light drives photochemistry
==================================================================

The first law of photochemistry (Grotthuss 1817, Draper 1842) says that
light the sample does not absorb produces no chemical change. The
quantity that matters is therefore the *absorbed* photon flux,
:math:`I_{abs}=I_0(1-10^{-A})`, computed by
:func:`~chemistrykit.photochem.photons_absorbed`. Here two samples with
the same incident flux but different absorbances form product in
proportion to the light each one actually absorbs, not the light that
falls on it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import photons_absorbed

I0 = 1.0e-6  # incident photon flux, einstein/s
absorbance = np.linspace(0.0, 3.0, 200)
I_abs = photons_absorbed(I0, absorbance)
I_trans = I0 - I_abs

fig, ax = plt.subplots()
ax.plot(absorbance, I_abs / I0, label="absorbed (acts)")
ax.plot(absorbance, I_trans / I0, "--", label="transmitted (does nothing)")
ax.set_xlabel("Absorbance $A$ at the irradiation wavelength")
ax.set_ylabel(r"Fraction of incident flux $I_0$")
ax.set_title("Grotthuss-Draper: only the absorbed fraction is photochemically active")
ax.legend()
fig.tight_layout()

# %%
# Same lamp, same exposure, same quantum yield: a weakly absorbing sample
# (A = 0.05) and a strongly absorbing one (A = 1.0). Product formed tracks
# the absorbed photons; a sample that absorbs nothing (A = 0) forms
# nothing, however bright the lamp.
Phi, t_exposure = 0.5, 600.0
for A in (0.0, 0.05, 1.0):
    n_product = Phi * photons_absorbed(I0, A) * t_exposure
    print(f"A = {A:4.2f}: absorbed fraction = {1 - 10**-A:.3f}, product formed = {n_product:.2e} mol")

plt.show()
