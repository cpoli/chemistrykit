r"""
Davy's electrolytic isolation of potassium and sodium
=======================================================

In 1807 Humphry Davy passed current through molten potash and soda and
isolated potassium and sodium -- metals no chemical reducing agent of the
day could free, because their reduction potentials are among the most
negative in the electrochemical series. This example ranks Davy's two
couples against the rest of
:data:`~chemistrykit.electrochem.systems.standard_potentials.STANDARD_REDUCTION_POTENTIALS`,
shows why the common reducing metals (Zn, Fe) cannot reduce K+ or Na+,
computes the minimum electrolysis voltage for the molten chlorides with
:func:`~chemistrykit.electrochem.systems.electrolysis.minimum_applied_voltage_electrolytic`,
and uses Faraday's law to estimate how much potassium a current yields.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.electrolysis import faradays_law_mass, minimum_applied_voltage_electrolytic
from chemistrykit.electrochem.systems.standard_potentials import STANDARD_REDUCTION_POTENTIALS as T
from chemistrykit.electrochem.systems.standard_potentials import cell_potential, is_spontaneous

# %%
# Could zinc or iron (the reducing agents available in 1807) reduce K+?
# The cell "K+ reduced by Zn" has a strongly negative potential: no.
for metal in ("Zn2+/Zn", "Fe2+/Fe"):
    for target in ("K+/K", "Na+/Na"):
        E = cell_potential(cathode=T[target], anode=T[metal])
        print(f"{target} reduced by {metal.split('/')[1]}: E = {E:+.2f} V, spontaneous = {is_spontaneous(E)}")

# %%
# Driving the reaction electrically: minimum voltage to decompose the
# chloride into metal (cathode) and chlorine (anode), using standard
# aqueous potentials as an estimate for the molten salt.
for target in ("K+/K", "Na+/Na"):
    E = cell_potential(cathode=T[target], anode=T["Cl2/Cl-"])
    print(f"Electrolysis of {target.split('/')[1]}Cl: at least {minimum_applied_voltage_electrolytic(E):.2f} V")

mass_K = faradays_law_mass(current=1.0, time=3600.0, molar_mass=39.098, n=1)
print(f"\n1 A for 1 h deposits {mass_K:.2f} g of potassium")

# %%
names = list(T)
values = np.array([T[name].E_standard for name in names])
order = np.argsort(values)
colors = ["crimson" if names[i] in ("K+/K", "Na+/Na") else "lightgray" for i in order]
fig, ax = plt.subplots(figsize=(6, 7))
ax.barh(range(len(names)), values[order], color=colors)
ax.set_yticks(range(len(names)))
ax.set_yticklabels([names[i] for i in order], fontsize=8)
ax.axvline(0.0, color="black", linewidth=0.8)
ax.set_xlabel(r"$E^\circ$ (V vs. SHE)")
ax.set_title("Davy's targets: K and Na near the bottom of the series")
fig.tight_layout()
plt.show()
