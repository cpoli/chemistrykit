r"""
Vavilov's law: fluorescence quantum yield is independent of excitation wavelength
===================================================================================

Vavilov found (1922-1927) that a dye's fluorescence quantum yield does not
change when it is excited at shorter, more energetic wavelengths. The
reason is that ultrafast internal conversion carries any higher excited
state down to :math:`S_1` before emission, so the yield depends only on the
:math:`S_1` rate constants,
:math:`\Phi_f=k_f/(k_f+k_{ic}+k_{isc})`
(:func:`~chemistrykit.photochem.fluorescence_quantum_yield`). Below, a
model dye with an :math:`S_0\to S_1` band near 480 nm and an
:math:`S_0\to S_2` band near 330 nm is excited across its whole
absorption spectrum; :func:`~chemistrykit.photochem.kasha_emission_yields`
gives the :math:`S_1` fluorescence yield for each state reached.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import fluorescence_quantum_yield, kasha_emission_yields

kf, kic, kisc = 2.0e8, 1.0e8, 0.5e8  # S1 rate constants, 1/s
kf2, k_ic21 = 1.0e8, 1.0e13  # S2 radiative and S2 -> S1 internal conversion
phi_S1 = fluorescence_quantum_yield(kf, kic, kisc)

wavelength = np.linspace(300.0, 520.0, 300)
band_S1 = np.exp(-(((wavelength - 480.0) / 20.0) ** 2))
band_S2 = 1.6 * np.exp(-(((wavelength - 330.0) / 18.0) ** 2))
frac_S2 = band_S2 / (band_S1 + band_S2)  # fraction of absorptions landing in S2

_, phi_via_S2 = kasha_emission_yields(kf2, k_ic21, kf, kic + kisc, excite="S2")
_, phi_via_S1 = kasha_emission_yields(kf2, k_ic21, kf, kic + kisc, excite="S1")
phi_f = frac_S2 * phi_via_S2 + (1.0 - frac_S2) * phi_via_S1

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))
ax1.plot(wavelength, band_S1 + band_S2)
ax1.set_ylabel("Absorbance (a.u.)")
ax1.set_title("Model dye absorption: $S_0\\to S_2$ (330 nm) and $S_0\\to S_1$ (480 nm)")
ax2.plot(wavelength, phi_f)
ax2.axhline(phi_S1, color="k", ls=":", label=r"$k_f/(k_f+k_{ic}+k_{isc})$")
ax2.set_ylim(0, 1)
ax2.set_xlabel("Excitation wavelength (nm)")
ax2.set_ylabel(r"$\Phi_f$")
ax2.set_title("Vavilov's law: flat fluorescence quantum yield")
ax2.legend()
fig.tight_layout()

print(f"Phi_f spans {phi_f.min():.6f} to {phi_f.max():.6f} across 300-520 nm (S1 value {phi_S1:.6f})")

plt.show()
