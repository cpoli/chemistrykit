r"""
Golay's open-tubular column: the van Deemter equation with A = 0
===================================================================

An open capillary has no packing, so Golay's equation for it is the van
Deemter equation with no eddy-diffusion term: :math:`H=B/u+Cu`. Using
:func:`~chemistrykit.analytical.van_deemter_H` with :math:`A=0` for a
capillary and :math:`A>0` for a packed column shows the lower plate
height, and hence the many more plates a long capillary delivers.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical import minimum_plate_height, optimum_flow_velocity, van_deemter_H

u = np.linspace(1.0, 100.0, 3000)  # cm/s
# Illustrative GC coefficients (H in cm): packed column vs. open-tubular capillary.
packed = {"A": 0.05, "B": 0.4, "C": 0.002, "length_cm": 200.0}
capillary = {"A": 0.0, "B": 0.4, "C": 0.0005, "length_cm": 3000.0}

fig, ax = plt.subplots(figsize=(7, 4.5))
for name, p in (("packed column (A > 0)", packed), ("open tubular, Golay (A = 0)", capillary)):
    H = van_deemter_H(u, p["A"], p["B"], p["C"])
    H_min = minimum_plate_height(p["A"], p["B"], p["C"])
    u_opt = optimum_flow_velocity(p["B"], p["C"])
    N = p["length_cm"] / H_min  # N = L / H
    print(f"{name:30s}: H_min = {H_min * 10:.3f} mm at u = {u_opt:.1f} cm/s; {p['length_cm'] / 100:.0f} m column -> N = {N:,.0f} plates")
    ax.plot(u, H * 10, label=name)
    ax.plot([u_opt], [H_min * 10], "o", color=ax.lines[-1].get_color())

# %%
ax.set_ylim(0, 2.0)
ax.set_xlabel("linear velocity u (cm/s)")
ax.set_ylabel("plate height H (mm)")
ax.set_title("Removing the eddy-diffusion term lowers H")
ax.legend()
plt.tight_layout()
plt.show()
