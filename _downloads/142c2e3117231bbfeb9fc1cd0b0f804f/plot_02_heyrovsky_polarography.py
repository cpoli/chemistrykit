r"""
Heyrovský's polarography and the Ilkovič equation
===================================================

Heyrovský's dropping mercury electrode (1922) records a sigmoidal
"polarographic wave" for each reducible species: its position (the
half-wave potential) identifies the species and its height (the
diffusion current, given by Ilkovič's 1934 equation) measures its
concentration. This example builds a two-component polarogram from
:func:`~chemistrykit.electrochem.systems.voltammetry.polarographic_wave_current`
with wave heights from
:func:`~chemistrykit.electrochem.systems.voltammetry.ilkovic_diffusion_current`,
and shows the linear calibration of diffusion current against concentration.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.voltammetry import ilkovic_diffusion_current, polarographic_wave_current

m, t_drop = 2.0, 4.0  # mercury flow (mg/s) and drop time (s)
# Two-electron reductions of Cd2+ and Zn2+ (D in cm^2/s, E1/2 in V vs. SCE, C in mmol/L).
species = {"Cd2+": (7.2e-6, -0.60, 0.5), "Zn2+": (7.0e-6, -1.00, 1.0)}

# %%
E = np.linspace(-0.3, -1.3, 400)
total = np.zeros_like(E)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
for name, (D, E_half, C) in species.items():
    i_d = ilkovic_diffusion_current(2, D, m, t_drop, C, average=True)
    total += polarographic_wave_current(E, E_half, i_d, n=2)
    print(f"{name}: E1/2 = {E_half:+.2f} V, C = {C} mM, i_d = {i_d:.2f} uA")
ax1.plot(E, total)
ax1.invert_xaxis()
ax1.set_xlabel("E (V vs. SCE)")
ax1.set_ylabel(r"Mean current ($\mu$A)")
ax1.set_title("Polarogram: one wave per reducible ion")

# %%
# Calibration: the Ilkovič diffusion current is proportional to concentration.
C_cal = np.linspace(0.1, 2.0, 8)
i_cal = [ilkovic_diffusion_current(2, species["Cd2+"][0], m, t_drop, c, average=True) for c in C_cal]
ax2.plot(C_cal, i_cal, "o-")
ax2.set_xlabel(r"[Cd$^{2+}$] (mmol/L)")
ax2.set_ylabel(r"$i_d$ ($\mu$A)")
ax2.set_title("Ilkovič calibration line")
fig.tight_layout()
plt.show()
