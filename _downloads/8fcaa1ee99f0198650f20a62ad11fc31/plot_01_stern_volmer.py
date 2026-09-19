r"""
Stern-Volmer quenching, and distinguishing static from dynamic mechanisms
=============================================================================

:func:`~chemistrykit.photochem.systems.stern_volmer.stern_volmer_ratio`
gives the linear intensity-ratio-vs-quencher-concentration relationship;
:func:`~chemistrykit.photochem.systems.stern_volmer.fit_stern_volmer`
recovers the Stern-Volmer constant from synthetic data. Measuring both
the intensity-ratio and lifetime-ratio slopes lets
:func:`~chemistrykit.photochem.systems.stern_volmer.classify_quenching_mechanism`
distinguish a purely dynamic (collisional) mechanism from a purely
static (ground-state complexation) one.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.photochem.systems.stern_volmer import (
    classify_quenching_mechanism,
    dynamic_quenching_constant,
    fit_stern_volmer,
    stern_volmer_ratio,
)
from chemistrykit.photochem.visualizers.photochem_plots import plot_stern_volmer

kq, tau0 = 2.0e10, 5.0e-9  # diffusion-controlled quenching, 5 ns unquenched lifetime
Ksv = dynamic_quenching_constant(kq, tau0)
print(f"Ksv (dynamic) = {Ksv:.3f} 1/M")

Q = np.array([0.0, 0.005, 0.010, 0.020, 0.040])
intensity_ratio = stern_volmer_ratio(Ksv, Q)
fit = fit_stern_volmer(Q, intensity_ratio)
print(f"Fitted Ksv = {fit.Ksv:.3f} 1/M, R^2 = {fit.r_squared:.6f}")

ax = plot_stern_volmer(Q, intensity_ratio, fit=fit)
plt.tight_layout()

# %%
# For purely dynamic quenching, the lifetime ratio tracks the intensity
# ratio exactly -- the two Stern-Volmer slopes agree.
lifetime_ratio = stern_volmer_ratio(Ksv, Q)  # dynamic: tau0/tau follows the same law as I0/I
mechanism_dynamic = classify_quenching_mechanism(intensity_ratio_slope=fit.Ksv, lifetime_ratio_slope=Ksv)
print(f"\nMechanism (equal slopes): {mechanism_dynamic}")

# For purely static quenching, the lifetime is unaffected (slope 0)
# even though the intensity ratio still rises linearly.
mechanism_static = classify_quenching_mechanism(intensity_ratio_slope=fit.Ksv, lifetime_ratio_slope=0.0)
print(f"Mechanism (zero lifetime slope): {mechanism_static}")

plt.show()
