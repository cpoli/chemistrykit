r"""
The Einstein vibrational heat capacity, and Sackur-Tetrode entropy
=====================================================================

:class:`~chemistrykit.statmech.VibrationalPartitionFunctionHarmonic`'s
heat capacity is the Einstein-oscillator formula: it vanishes as
:math:`T\to0` (the vibration is "frozen out") and approaches the
classical equipartition value :math:`k_B` per mode as
:math:`T\to\infty`. Separately,
:func:`~chemistrykit.statmech.sackur_tetrode_entropy` gives the
translational entropy of an ideal gas directly from its molecular mass,
temperature, and pressure -- checked here against the standard tabulated
value for argon.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

from chemistrykit.statmech import VibrationalPartitionFunctionHarmonic, sackur_tetrode_entropy
from chemistrykit.statmech.visualizers.statmech_plots import plot_heat_capacity_vs_temperature

K_B = sc.k

hcl_mode = VibrationalPartitionFunctionHarmonic.from_wavenumber(2886.0)  # HCl fundamental, cm^-1
print(f"HCl vibrational temperature: {hcl_mode.vibrational_temperature:.1f} K")

T = np.linspace(50.0, 6000.0, 300)

fig, ax = plt.subplots(figsize=(7, 5))
plot_heat_capacity_vs_temperature(hcl_mode, T, ax=ax, color="steelblue")
ax.axhline(K_B, color="gray", linestyle=":", linewidth=0.8, label="classical limit (k_B)")
ax.axvline(hcl_mode.vibrational_temperature, color="crimson", linestyle=":", linewidth=0.8, label="vibrational temperature")
ax.set_title("Einstein vibrational heat capacity (HCl-like mode)")
ax.legend()
fig.tight_layout()

# %%
# At room temperature, HCl's vibrational mode is still almost entirely
# frozen out (T << vibrational temperature), so it contributes almost
# nothing to the heat capacity -- the classical equipartition
# expectation only becomes accurate at temperatures far higher than any
# real HCl gas would survive at:

Cv_room = hcl_mode.heat_capacity_v(298.15, N=1.0) / K_B
print(f"Cv_vib(298.15 K) / k_B = {Cv_room:.4f} (classical limit would be 1.0)")

# %%
# Separately, the Sackur-Tetrode equation gives argon's standard molar
# translational entropy directly from its mass -- matching the tabulated
# value (Atkins & de Paula, Table 13.1) to three figures:

argon_mass = 39.948 * sc.atomic_mass
S = sackur_tetrode_entropy(mass=argon_mass, T=298.15, P=1.0e5)
print(f"S_m(Ar, 298.15 K, 1 bar) = {S:.2f} J/(mol K)  (tabulated: ~154.8 J/(mol K))")

plt.show()
