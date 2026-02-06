#!/usr/bin/env python3
"""
Generate Table S2: Pangenome ortholog group counts with/without paralog splitting.
Reads Panaroo summary_statistics.txt files and performs proportions z-tests
and a chi-square test.

Usage:
    python generate_table_s2.py <split_paralog_file> <no_split_paralog_file>

Example:
    python generate_table_s2.py \
        No_merge_paralogs_summary_statistics.txt \
        merge_paralogs_summary_statistics.txt
"""

import sys
import numpy as np
from scipy import stats


def parse_summary(filepath):
    """Parse a Panaroo summary_statistics.txt file into a dict of counts."""
    mapping = {
        "Core genes": "core",
        "Soft core genes": "soft_core",
        "Shell genes": "shell",
        "Cloud genes": "cloud",
        "Total genes": "total",
    }
    counts = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3:
                label = parts[0].strip()
                count = int(parts[2].strip())
                if label in mapping:
                    counts[mapping[label]] = count
    return counts


def proportions_ztest(x1, n1, x2, n2):
    """Two-sided z-test for equality of two proportions."""
    p1 = x1 / n1
    p2 = x2 / n2
    p_pool = (x1 + x2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / se
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p_val


def significance_stars(p):
    """Return significance stars for a p-value."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    split_file = sys.argv[1]      # No_merge_paralogs (split)
    no_split_file = sys.argv[2]   # merge_paralogs (no split)

    split = parse_summary(split_file)
    no_split = parse_summary(no_split_file)

    categories = [
        ("Core (99%+)", "core"),
        ("Soft core (95–99%)", "soft_core"),
        ("Shell (15–95%)", "shell"),
        ("Cloud (<15%)", "cloud"),
    ]

    n_split = split["total"]
    n_nosplit = no_split["total"]

    # Header
    header = (
        f"{'Genes (% of isolates)':<25} "
        f"{'Split paralog':>15} "
        f"{'No split paralog':>18} "
        f"{'Sig.':>5} "
        f"{'P value':>12} "
        f"{'Test':<20}"
    )
    sep = "-" * len(header)

    print(sep)
    print(header)
    print(sep)

    # Per-category z-tests
    for label, key in categories:
        z, p = proportions_ztest(split[key], n_split, no_split[key], n_nosplit)
        stars = significance_stars(p)
        print(
            f"{label:<25} "
            f"{split[key]:>15} "
            f"{no_split[key]:>18} "
            f"{stars:>5} "
            f"{p:>12.2e} "
            f"{'proportions z-test':<20}"
        )

    # Chi-square on full 2x4 table
    table = np.array([
        [split[k] for _, k in categories],
        [no_split[k] for _, k in categories],
    ])
    chi2, p_chi, dof, _ = stats.chi2_contingency(table)
    stars = significance_stars(p_chi)
    print(
        f"{'Total':<25} "
        f"{n_split:>15} "
        f"{n_nosplit:>18} "
        f"{stars:>5} "
        f"{p_chi:>12.2e} "
        f"{f'χ² test (df={dof})':<20}"
    )
    print(sep)

    print(
        "\nTable S2: Counts of ortholog groups in the pangenome, with and "
        "without splitting of paralogous gene families. The p value reports "
        "the test of hypothesis that counts of split and unsplit paralogs "
        "are the same, using a z-test of proportions."
    )


if __name__ == "__main__":
    main()
