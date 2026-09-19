r"""
Dixon's Q-test for outlier rejection
========================================

:func:`~chemistrykit.analytical.systems.qtest.dixon_q_test` compares the
gap between a suspect value and its nearest neighbor, relative to the
full data range, against a tabulated critical value
(:data:`~chemistrykit.analytical.systems.qtest.Q_CRITICAL_TABLE`) that
depends on sample size and confidence level.
"""

# %%
from chemistrykit.analytical.systems.qtest import Q_CRITICAL_TABLE, dixon_q_test

# %%
# Five replicate titration endpoints (mL), one of them clearly off:
replicates = [24.51, 24.55, 24.48, 24.53, 25.10]

for confidence in (0.90, 0.95, 0.99):
    result = dixon_q_test(replicates, confidence=confidence)
    verdict = "REJECT" if result.reject else "retain"
    print(f"confidence={confidence:.0%}: suspect={result.suspect_value}, Q={result.Q_statistic:.4f}, Q_crit={result.Q_critical:.4f} -> {verdict}")

# %%
# The critical value grows with confidence level (harder to reject at
# higher confidence) and shrinks with sample size (more data makes a
# single outlier easier to identify):
print("\nQ_crit(n, 95%) for n=3..10:")
for n in range(3, 11):
    print(f"  n={n}: {Q_CRITICAL_TABLE[n][0.95]}")

# %%
# If the suspect value were retained instead (a smaller, more plausible
# deviation), the Q-test would not reject it even at 90% confidence:
close_replicates = [24.51, 24.55, 24.48, 24.53, 24.62]
result = dixon_q_test(close_replicates, confidence=0.90)
print(f"\nCloser data set: Q={result.Q_statistic:.4f}, Q_crit={result.Q_critical:.4f}, reject={result.reject}")
