r"""
Norrish and Porter's flash photolysis: watching a transient triplet decay
===========================================================================

Flash photolysis (Norrish and Porter, 1949) excites a sample with an
intense, short flash and then records its absorption at a series of
delays with a second probe flash, so a short-lived intermediate is seen
appearing and disappearing in real time. Here the intermediate is the
triplet state of the Jablonski model
(:func:`~chemistrykit.photochem.jablonski_populations_analytic`): the
probe measures the triplet-triplet absorbance
:math:`\Delta A(t)=\varepsilon_T\ell[T_1](t)`, and fitting the late-time
decay of noisy "measured" points recovers the triplet lifetime.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem import jablonski_populations_analytic

kf, kic, kisc = 1.0e8, 5.0e7, 3.5e8  # S1 decays in a few ns
kp, kic_T = 5.0e2, 1.95e4  # T1 lifetime = 50 microseconds
conc_excited = 2.0e-5  # M of molecules excited by the flash
eps_T, path = 1.0e4, 1.0  # triplet-triplet molar absorptivity (1/(M cm)), path (cm)

t_fine = np.linspace(0.0, 300e-6, 2000)
_, T1, _ = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, conc_excited, t_fine)
dA_true = eps_T * path * T1

rng = np.random.default_rng(0)
delays = np.linspace(10e-6, 250e-6, 25)  # probe-flash delays
_, T1_d, _ = jablonski_populations_analytic(kf, kic, kisc, kp, kic_T, conc_excited, delays)
dA_meas = eps_T * path * T1_d + rng.normal(0.0, 0.002, delays.size)

# %%
# Fit ln(dA) vs delay: slope = -1/tau_T (S1 is long gone after 10 us).
mask = dA_meas > 0
slope, intercept = np.polyfit(delays[mask], np.log(dA_meas[mask]), 1)
tau_fit = -1.0 / slope
print(f"True triplet lifetime: {1 / (kp + kic_T) * 1e6:.1f} us, fitted from probe delays: {tau_fit * 1e6:.1f} us")

fig, ax = plt.subplots()
ax.plot(t_fine * 1e6, dA_true, label=r"triplet absorbance $\varepsilon_T\ell[T_1](t)$")
ax.plot(delays * 1e6, dA_meas, "o", label="probe flash at each delay")
ax.plot(t_fine * 1e6, np.exp(intercept + slope * t_fine), "k--", label=rf"fit, $\tau_T$ = {tau_fit * 1e6:.0f} $\mu$s")
ax.set_xlabel(r"Delay after photolysis flash ($\mu$s)")
ax.set_ylabel(r"Transient absorbance $\Delta A$")
ax.set_title("Flash photolysis of a transient triplet")
ax.legend()
fig.tight_layout()

plt.show()
