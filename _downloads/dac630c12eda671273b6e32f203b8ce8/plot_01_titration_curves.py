r"""
Acid-base, redox, and EDTA titration curves side by side
============================================================

Three different titration *types* -- acid-base (pH,
:mod:`chemistrykit.solutions.systems.titration`), redox (electrode
potential `E`, :class:`~chemistrykit.analytical.systems.titration.RedoxTitration`),
and complexometric (`pM`,
:class:`~chemistrykit.analytical.systems.titration.EDTATitration`) --
all share the same qualitative S-shaped curve, and the same numerical
equivalence-point detector (steepest response change).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.analytical.systems.titration import EDTATitration, RedoxTitration
from chemistrykit.analytical.visualizers.analytical_plots import plot_titration_curve
from chemistrykit.solutions.systems.titration import WeakAcidStrongBaseTitration

# %%
# Acid-base: acetic acid titrated with NaOH (reusing chemistrykit.solutions).
acid_base = WeakAcidStrongBaseTitration(Ca=0.100, Va=0.050, Ka=1.8e-5, Cb=0.100)
Vb_ab = np.linspace(1e-6, 0.099, 4000)
curve_ab = acid_base.curve(Vb_ab)
Veq_ab_numeric = acid_base.find_equivalence_point(Vb_ab)
print(f"Acid-base equivalence volume: exact={acid_base.equivalence_volume() * 1000:.3f} mL, numeric={Veq_ab_numeric * 1000:.3f} mL")

# %%
# Redox: Fe2+ titrated with Ce4+ (both 1-electron couples).
redox = RedoxTitration(E1_standard=0.771, n1=1, E2_standard=1.72, n2=1, C_analyte=0.100, V_analyte=0.050, C_titrant=0.100)
V_redox = np.linspace(1e-6, 0.099, 4000)
curve_redox = redox.curve(V_redox)
Veq_redox_numeric = redox.find_equivalence_point(V_redox)
print(f"Redox equivalence volume:     exact={redox.equivalence_volume() * 1000:.3f} mL, numeric={Veq_redox_numeric * 1000:.3f} mL")
print(f"E at equivalence: {float(redox.response_at(np.array([redox.equivalence_volume()]))[0]):.4f} V (= arithmetic mean of the two standard potentials)")

# %%
# Complexometric: Ca2+ titrated with EDTA at a fixed pH (conditional
# formation constant K'=10^10, typical for EDTA-Ca2+ around pH 10).
edta = EDTATitration(C_metal=0.0100, V_metal=0.050, K_conditional=1e10, C_edta=0.0100)
V_edta = np.linspace(1e-6, 0.099, 4000)
curve_edta = edta.curve(V_edta)
Veq_edta_numeric = edta.find_equivalence_point(V_edta)
print(f"EDTA equivalence volume:      exact={edta.equivalence_volume() * 1000:.3f} mL, numeric={Veq_edta_numeric * 1000:.3f} mL")

# %%
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
plot_titration_curve(curve_ab.Vb * 1000, curve_ab.pH, V_equiv=acid_base.equivalence_volume() * 1000, ax=axes[0], ylabel="pH")
axes[0].set_title("Acid-base (weak acid + NaOH)")
plot_titration_curve(curve_redox.V * 1000, curve_redox.response, V_equiv=redox.equivalence_volume() * 1000, ax=axes[1], ylabel="E (V)")
axes[1].set_title("Redox (Fe2+ + Ce4+)")
plot_titration_curve(curve_edta.V * 1000, curve_edta.response, V_equiv=edta.equivalence_volume() * 1000, ax=axes[2], ylabel="pCa")
axes[2].set_title("Complexometric (Ca2+ + EDTA)")
for ax in axes:
    ax.set_xlabel("titrant volume (mL)")
plt.tight_layout()
plt.show()
