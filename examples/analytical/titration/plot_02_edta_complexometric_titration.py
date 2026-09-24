r"""
EDTA complexometric titration of a metal ion
===============================================

:class:`~chemistrykit.analytical.EDTATitration` follows
:math:`M+Y\rightleftharpoons MY` through its conditional formation
constant :math:`K_f'`. The same EDTA titrant quantifies any metal it
binds strongly enough; the size of the `pM` break at equivalence is set
by :math:`K_f'`, which is why Schwarzenbach's titrations are run in a
buffer at a pH where :math:`K_f'` is large.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical import EDTATitration

V = np.linspace(1e-6, 0.080, 4000)
curves = {}
for logK in (6, 8, 10, 12):
    titration = EDTATitration(C_metal=0.0100, V_metal=0.050, K_conditional=10.0**logK, C_edta=0.0100)
    curves[logK] = titration.curve(V)
    V_eq = titration.equivalence_volume()
    pM_before, pM_eq, pM_after = titration.response_at(np.array([0.99 * V_eq, V_eq, 1.01 * V_eq]))
    C_M_eq = 0.0100 * 0.050 / (0.050 + V_eq)
    print(
        f"log K'={logK:2d}: pM at equivalence {pM_eq:.3f} (large-K approximation {0.5 * np.log10(10.0**logK / C_M_eq):.3f}), "
        f"jump over +/-1% of V_eq = {pM_after - pM_before:.2f} pM units, endpoint found at {titration.find_equivalence_point(V) * 1000:.2f} mL"
    )

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
for logK, curve in curves.items():
    ax.plot(curve.V * 1000, curve.response, label=f"log $K_f'$ = {logK}")
ax.axvline(50.0, color="gray", linestyle="--", linewidth=0.8)
ax.set_xlabel("EDTA added (mL)")
ax.set_ylabel("pM")
ax.set_title("Complexometric titration of 0.0100 M M$^{2+}$ with 0.0100 M EDTA")
ax.legend()
plt.tight_layout()
plt.show()
