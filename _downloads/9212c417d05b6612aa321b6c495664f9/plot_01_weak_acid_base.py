r"""
The Ostwald dilution law and buffer design
============================================

Dissolving less acid in more water increases the *percentage* of it that
dissociates, even though the absolute :math:`[H^+]` goes down -- the
Ostwald dilution law, here from the exact cubic charge-balance solution
in :class:`~chemistrykit.solutions.systems.acid_base.WeakAcid` rather
than the more commonly taught :math:`\sqrt{K_aC_a}` approximation. Also:
designing a buffer for a target pH with
:meth:`~chemistrykit.solutions.systems.acid_base.Buffer.from_target_ph`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.solutions.systems.acid_base import Buffer, WeakAcid

Ka = 1.8e-5  # acetic acid
Ca_values = np.logspace(-4, 0, 60)
percent_dissociated = [WeakAcid(Ca=Ca, Ka=Ka).percent_dissociation() for Ca in Ca_values]

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogx(Ca_values, percent_dissociated, color="steelblue")
ax.set_xlabel("total acid concentration Ca (mol/L)")
ax.set_ylabel("percent dissociated")
ax.set_title("Ostwald dilution law: acetic acid")
fig.tight_layout()

# %%
# Designing an acetate buffer (pKa = 4.76) at pH 5.0 from a fixed total
# buffer concentration:

buf = Buffer.from_target_ph(pKa=4.76, target_pH=5.0, total_conc=0.20)
print(f"Buffer: [HA] = {buf.acid_conc:.4f} M, [A-] = {buf.base_conc:.4f} M, pH = {buf.pH():.4f}")

pH_targets = np.linspace(3.76, 5.76, 100)
base_fraction = [Buffer.from_target_ph(pKa=4.76, target_pH=pH, total_conc=0.20).base_conc / 0.20 for pH in pH_targets]

fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.plot(pH_targets, base_fraction, color="darkorange")
ax2.axvline(4.76, color="gray", linestyle=":", label="pH = pKa")
ax2.set_xlabel("target buffer pH")
ax2.set_ylabel("fraction as conjugate base, [A-]/C")
ax2.set_title("Buffer composition vs. target pH")
ax2.legend()
fig2.tight_layout()

plt.show()
