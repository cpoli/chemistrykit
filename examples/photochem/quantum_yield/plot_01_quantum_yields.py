r"""
Fluorescence/phosphorescence quantum yields, and photon absorption via Beer-Lambert
========================================================================================

:func:`~chemistrykit.photochem.systems.quantum_yield.fluorescence_quantum_yield`
and :func:`~chemistrykit.photochem.systems.quantum_yield.phosphorescence_quantum_yield`
are simple branching-ratio calculations from Jablonski-diagram rate
constants. :func:`~chemistrykit.photochem.systems.quantum_yield.photons_absorbed`
reuses :mod:`chemistrykit.spectro.systems.beer_lambert`'s transmittance
to compute the photon flux actually absorbed by a sample, which feeds
directly into :func:`~chemistrykit.photochem.systems.quantum_yield.photochemical_quantum_yield`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem.systems.quantum_yield import (
    fluorescence_quantum_yield,
    phosphorescence_quantum_yield,
    photochemical_quantum_yield,
    photons_absorbed,
)

# %%
# Fluorescence quantum yield as internal conversion competes more and
# more strongly against fluorescence.
kf = 1.0
kic_values = np.linspace(0.0, 5.0, 50)
kisc = 0.2
phi_f = np.array([fluorescence_quantum_yield(kf, kic, kisc) for kic in kic_values])

fig, ax = plt.subplots()
ax.plot(kic_values, phi_f)
ax.set_xlabel(r"$k_{ic}$")
ax.set_ylabel(r"$\Phi_f$")
ax.set_title("Fluorescence quantum yield vs. competing internal conversion")
fig.tight_layout()

# %%
# Phosphorescence quantum yield requires *both* successful intersystem
# crossing *and* a radiative triplet decay -- the product of two
# branching ratios, so it is always <= either one alone.
kf, kic, kp, kic_T = 1.0, 1.0, 0.5, 0.5
kisc_values = np.linspace(0.0, 3.0, 50)
phi_p = np.array([phosphorescence_quantum_yield(kisc, kf, kic, kp, kic_T) for kisc in kisc_values])

fig2, ax2 = plt.subplots()
ax2.plot(kisc_values, phi_p)
ax2.set_xlabel(r"$k_{isc}$")
ax2.set_ylabel(r"$\Phi_p$")
ax2.set_title("Phosphorescence quantum yield vs. intersystem crossing rate")
fig2.tight_layout()

# %%
# Photochemical quantum yield: photons absorbed via Beer-Lambert, then
# product formed per photon absorbed.
I0 = 1.0e-6  # mol photons / s
absorbance = 0.8
absorbed = photons_absorbed(I0, absorbance)
moles_product_per_second = 0.4 * absorbed  # a hypothetical measured product-formation rate
Phi = photochemical_quantum_yield(moles_product_per_second, absorbed)
print(f"Photon flux absorbed: {absorbed:.3e} mol/s")
print(f"Photochemical quantum yield: {Phi:.3f}")

plt.show()
