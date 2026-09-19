r"""
Plate theory and the van Deemter equation
=============================================

:func:`~chemistrykit.analytical.systems.chromatography.theoretical_plates`
converts a peak's retention time and width into a plate count;
:func:`~chemistrykit.analytical.systems.chromatography.van_deemter_H`
then predicts how plate height varies with flow velocity, with a minimum
(the optimum operating velocity) at :math:`u_{opt}=\sqrt{B/C}`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical.systems.chromatography import (
    minimum_plate_height,
    optimum_flow_velocity,
    plate_height,
    resolution,
    selectivity_factor,
    simulate_chromatogram,
    theoretical_plates,
    van_deemter_H,
)
from chemistrykit.analytical.visualizers.analytical_plots import plot_chromatogram, plot_van_deemter

# %%
# A peak eluting at 10.0 min with a 0.4 min base width:
N = theoretical_plates(retention_time=10.0, peak_width=0.4, width_type="base")
H = plate_height(column_length=25.0, N=N)
print(f"N = {N:.0f} theoretical plates")
print(f"H = {H * 1e4:.2f} um (25 cm column)")

# %%
# The van Deemter equation predicts an optimum flow velocity that
# minimizes plate height -- verified against a fine numerical scan of
# dH/du:
A, B, C = 1.5, 25.0, 0.05
u_opt = optimum_flow_velocity(B, C)
H_min = minimum_plate_height(A, B, C)
u = np.linspace(0.5, 40.0, 2000)
H_curve = van_deemter_H(u, A, B, C)
u_numeric = u[np.argmin(H_curve)]
print(f"\nu_opt (closed form) = {u_opt:.4f}, numeric scan minimum at u = {u_numeric:.4f}")
print(f"H_min (closed form) = {H_min:.4f}, numeric scan minimum H = {H_curve.min():.4f}")

# %%
# Resolution and selectivity between two closely-eluting peaks:
tR1, tR2 = 9.6, 10.4
w1, w2 = 0.35, 0.38
Rs = resolution(tR1, tR2, w1, w2)
alpha = selectivity_factor(k1=(tR1 - 1.0) / 1.0, k2=(tR2 - 1.0) / 1.0)
print(f"\nResolution Rs = {Rs:.3f} ({'baseline-separated' if Rs >= 1.5 else 'not fully resolved'})")
print(f"Selectivity alpha = {alpha:.4f}")

# %%
t = np.linspace(8.0, 12.0, 4000)
chrom = simulate_chromatogram(t, centers=[tR1, tR2], N=N)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
plot_van_deemter(u, H_curve, u_opt=u_opt, H_min=H_min, ax=axes[0])
plot_chromatogram(t, chrom, ax=axes[1])
plt.tight_layout()
plt.show()
