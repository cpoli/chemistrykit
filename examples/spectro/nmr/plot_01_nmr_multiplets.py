r"""
First-order NMR multiplets: ethanol's triplet/quartet, and a doublet of triplets
====================================================================================

Simulates ethanol's classic -CH3 triplet and -CH2- quartet (each proton
set split only by its 3-bond neighbor, the textbook first-order pattern),
then a genuine doublet of triplets from two chemically distinct coupling
partners with different `J` values.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.spectro.systems.nmr import first_order_multiplet, multi_coupling_multiplet
from chemistrykit.spectro.visualizers.spectro_plots import plot_broadened_spectrum, plot_stick_spectrum

freq_mhz = 400.0
ch3 = first_order_multiplet(chemical_shift_ppm=1.2, j_coupling_hz=7.0, n_neighbors=2, spectrometer_frequency_mhz=freq_mhz)  # split by CH2's 2 H's -> triplet
ch2 = first_order_multiplet(chemical_shift_ppm=3.7, j_coupling_hz=7.0, n_neighbors=3, spectrometer_frequency_mhz=freq_mhz)  # split by CH3's 3 H's -> quartet

print(f"CH3 triplet positions (ppm): {np.round(ch3.positions, 4)}, intensities: {ch3.intensities}")
print(f"CH2 quartet positions (ppm): {np.round(ch2.positions, 4)}, intensities: {ch2.intensities}")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
plot_stick_spectrum(ch3, ax=axes[0], color="steelblue")
axes[0].set_title("Ethanol -CH3 (triplet)")
axes[0].set_xlabel("chemical shift (ppm)")
plot_stick_spectrum(ch2, ax=axes[1], color="crimson")
axes[1].set_title("Ethanol -CH2- (quartet)")
axes[1].set_xlabel("chemical shift (ppm)")
fig.tight_layout()

# %%
# A genuine doublet of triplets: one neighbor with J1=12 Hz, two
# equivalent neighbors with J2=5 Hz -- 6 distinct lines, total relative
# intensity conserved at (1+1)*(1+2+1)=8:

dt = multi_coupling_multiplet(chemical_shift_ppm=5.5, couplings=[(12.0, 1), (5.0, 2)], spectrometer_frequency_mhz=freq_mhz)
print(f"\nDoublet-of-triplets: {len(dt.positions)} lines, total intensity {np.sum(dt.intensities):.1f}")
print(f"Positions (ppm): {np.round(dt.positions, 4)}")
print(f"Intensities: {dt.intensities}")

fig2, ax2 = plt.subplots(figsize=(7, 4))
x = np.linspace(dt.positions[0] - 0.02, dt.positions[-1] + 0.02, 3000)
plot_broadened_spectrum(dt, x, ax=ax2, shape="lorentzian", fwhm=0.002, color="darkgreen")
ax2.set_xlabel("chemical shift (ppm)")
ax2.set_title("Doublet of triplets (J1=12 Hz, J2=5 Hz)")
fig2.tight_layout()
plt.show()
