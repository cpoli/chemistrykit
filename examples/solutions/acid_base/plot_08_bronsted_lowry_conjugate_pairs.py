r"""
Bronsted-Lowry conjugate acid-base pairs: pKa + pKb = pKw
===========================================================

In the Bronsted-Lowry picture an acid is a proton donor and its
deprotonated form is a *conjugate base* -- itself a genuine base that
accepts a proton from water. For every pair, :math:`K_aK_b = K_w`, so the
stronger the acid, the weaker its conjugate base. Here the acetate ion is
modelled as a base in its own right with
:class:`~chemistrykit.solutions.systems.acid_base.WeakBase`
(:math:`K_b = K_w/K_a`), alongside acetic acid as
:class:`~chemistrykit.solutions.systems.acid_base.WeakAcid`, and several
conjugate pairs are placed on the line :math:`pK_a + pK_b = 14`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import WeakAcid, WeakBase

Kw = 1.0e-14
pairs = {
    "HF / F$^-$": 6.8e-4,
    "CH$_3$COOH / CH$_3$COO$^-$": 1.8e-5,
    "H$_2$CO$_3$ / HCO$_3^-$": 4.5e-7,
    "NH$_4^+$ / NH$_3$": 5.6e-10,
    "HCN / CN$^-$": 6.2e-10,
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
pKa_line = np.linspace(0, 14, 2)
ax1.plot(pKa_line, 14 - pKa_line, color="gray", linestyle="--", label="pKa + pKb = 14")
for name, Ka in pairs.items():
    pKa, pKb = -np.log10(Ka), -np.log10(Kw / Ka)
    ax1.plot(pKa, pKb, "o")
    ax1.annotate(name, (pKa, pKb), textcoords="offset points", xytext=(6, 2), fontsize=8)
ax1.set_xlabel("pKa of the acid")
ax1.set_ylabel("pKb of its conjugate base")
ax1.set_title("Conjugate pairs: stronger acid, weaker base")
ax1.legend()

# 0.10 M solutions of each acid and of its conjugate base (e.g. as the
# sodium salt): the acid solution is acidic, the conjugate-base solution
# is basic, and both are computed with the same exact cubic.

names = list(pairs)
pH_acid = [WeakAcid(Ca=0.10, Ka=Ka).pH() for Ka in pairs.values()]
pH_base = [WeakBase(Cb=0.10, Kb=Kw / Ka).pH() for Ka in pairs.values()]
x = np.arange(len(names))
ax2.bar(x - 0.2, pH_acid, width=0.4, color="indianred", label="0.10 M acid")
ax2.bar(x + 0.2, pH_base, width=0.4, color="steelblue", label="0.10 M conjugate base")
ax2.axhline(7.0, color="gray", linestyle=":")
ax2.set_xticks(x)
ax2.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
ax2.set_ylabel("pH")
ax2.set_title("pH of each member of the pair")
ax2.legend(fontsize=8)
fig.tight_layout()

# %%
# Acetate as a Bronsted base in its own right:

acetate = WeakBase(Cb=0.10, Kb=Kw / 1.8e-5)
print(f"0.10 M sodium acetate: Kb = {acetate.Kb:.3e}, pH = {acetate.pH():.3f}")

plt.show()
