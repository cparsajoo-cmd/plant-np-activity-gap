"""Exact power of Fisher's exact test for the prospective validation design.

This is an exact enumeration, not a simulation: every possible pair of outcome
counts is weighted by its binomial probability and tested, so the result is
deterministic and anyone running it obtains the identical number. No random
seed is involved.

Design: n1 compounds drawn from the model's prioritised set, n2 control
compounds, confirmation probabilities p1 and p2, Fisher's exact test, two
sided, alpha = 0.05.

Usage:  python3 04_power.py
"""
from math import comb
from scipy.stats import fisher_exact

ALPHA = 0.05


def binom_pmf(n, k, p):
    return comb(n, k) * p ** k * (1.0 - p) ** (n - k)


def power(n1, n2, p1, p2, alpha=ALPHA):
    """Exact probability that Fisher's exact test rejects, under p1 vs p2."""
    # Pre-compute which tables are significant; the test depends only on counts.
    total = 0.0
    w1 = [binom_pmf(n1, x, p1) for x in range(n1 + 1)]
    w2 = [binom_pmf(n2, y, p2) for y in range(n2 + 1)]
    for x in range(n1 + 1):
        if w1[x] == 0.0:
            continue
        for y in range(n2 + 1):
            if w2[y] == 0.0:
                continue
            p = fisher_exact([[x, n1 - x], [y, n2 - y]], alternative="two-sided")[1]
            if p <= alpha:
                total += w1[x] * w2[y]
    return total


DESIGNS = [(80, 80), (60, 60), (60, 30), (100, 100), (40, 40)]
EFFECTS = [(0.25, 0.05), (0.20, 0.05), (0.15, 0.05), (0.30, 0.05)]

print("Exact power of Fisher's exact test, two-sided, alpha = %.2f" % ALPHA)
print("(exhaustive enumeration over all outcome tables; no simulation)")
print()
header = "%-14s" % "design" + "".join("  %5.0f%% vs %2.0f%%" % (a * 100, b * 100) for a, b in EFFECTS)
print(header)
print("-" * len(header))
for n1, n2 in DESIGNS:
    row = "%-14s" % ("n=%d vs %d" % (n1, n2))
    for p1, p2 in EFFECTS:
        row += "  %12.1f%%" % (100.0 * power(n1, n2, p1, p2))
    print(row, flush=True)
print()
print("The design adopted in the proposal is n = 80 vs 80 (160 compounds in")
print("total), read from the first row.")
