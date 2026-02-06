#!/usr/bin/env python3
"""
Compare two Panaroo output directories by gene member overlap.

Matches genes between runs based on shared locus_tags (member overlap),
since group names are arbitrary and differ between runs.

Usage:
  uv run python3 panaroo/scripts/compare_panaroo_runs.py <dir1> <dir2> [--names name1 name2]

Example:
  uv run python3 panaroo/scripts/compare_panaroo_runs.py \
    panaroo/debugging_2/no_merge_paralogs \
    panaroo/output/panaroo_out \
    --names "my_run" "authors_run"
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict


def load_summary_stats(panaroo_dir: Path) -> dict:
    """Load summary_statistics.txt into a dict."""
    stats_file = panaroo_dir / "summary_statistics.txt"
    stats = {}
    with open(stats_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                category = parts[0].split()[0]
                count = int(parts[-1])
                stats[category] = count
    return stats


def load_groups_to_members(panaroo_dir: Path) -> dict:
    """Load gene_presence_absence.csv and extract group -> member locus_tags."""
    gpa_file = panaroo_dir / "gene_presence_absence.csv"
    df = pd.read_csv(gpa_file)

    sample_cols = df.columns[3:]  # Skip Gene, Non-unique Gene name, Annotation

    groups = {}
    for _, row in df.iterrows():
        group_id = row['Gene']
        annotation = row.get('Annotation', '')
        members = set()
        for col in sample_cols:
            val = row[col]
            if pd.notna(val) and val != '':
                genes = str(val).split(';')
                members.update([g.strip() for g in genes if g.strip()])
        groups[group_id] = {'members': members, 'annotation': annotation}

    return groups


def compare_stats(stats1: dict, stats2: dict, name1: str, name2: str):
    """Compare summary statistics."""
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS COMPARISON")
    print("=" * 60)

    print(f"\n{'Category':<15} {name1:>15} {name2:>15} {'Diff':>10}")
    print("-" * 55)

    for cat in ['Core', 'Soft', 'Shell', 'Cloud', 'Total']:
        v1 = stats1.get(cat, 0)
        v2 = stats2.get(cat, 0)
        diff = v1 - v2
        diff_str = f"+{diff}" if diff > 0 else str(diff)
        print(f"{cat:<15} {v1:>15} {v2:>15} {diff_str:>10}")


def match_groups_by_members(groups1: dict, groups2: dict, name1: str, name2: str):
    """Match groups between runs based on member overlap."""
    print("\n" + "=" * 60)
    print("GENE MATCHING BY MEMBER OVERLAP")
    print("=" * 60)

    # Build reverse index: locus_tag -> group for run2
    tag_to_group2 = {}
    for group, data in groups2.items():
        for tag in data['members']:
            tag_to_group2[tag] = group

    # Match groups from run1 to run2
    matches = {}  # group1 -> (group2, overlap_ratio)
    unmatched1 = []

    for group1, data1 in groups1.items():
        members1 = data1['members']

        # Find which groups in run2 share members
        group2_overlaps = defaultdict(int)
        for tag in members1:
            if tag in tag_to_group2:
                group2_overlaps[tag_to_group2[tag]] += 1

        if group2_overlaps:
            # Best match is the group with most shared members
            best_group2 = max(group2_overlaps, key=group2_overlaps.get)
            shared = group2_overlaps[best_group2]
            members2 = groups2[best_group2]['members']

            # Jaccard similarity
            union = len(members1 | members2)
            jaccard = shared / union if union else 0

            matches[group1] = {
                'group2': best_group2,
                'shared': shared,
                'only1': len(members1 - members2),
                'only2': len(members2 - members1),
                'jaccard': jaccard,
                'ann1': data1['annotation'],
                'ann2': groups2[best_group2]['annotation']
            }
        else:
            unmatched1.append(group1)

    # Find unmatched in run2
    matched_groups2 = set(m['group2'] for m in matches.values())
    unmatched2 = [g for g in groups2 if g not in matched_groups2]

    # Categorize matches
    perfect = [(g1, m) for g1, m in matches.items() if m['jaccard'] == 1.0]
    good = [(g1, m) for g1, m in matches.items() if 0.8 <= m['jaccard'] < 1.0]
    partial = [(g1, m) for g1, m in matches.items() if 0.5 <= m['jaccard'] < 0.8]
    poor = [(g1, m) for g1, m in matches.items() if m['jaccard'] < 0.5]

    print(f"\n{name1}: {len(groups1)} groups")
    print(f"{name2}: {len(groups2)} groups")

    print(f"\nMatch quality:")
    print(f"  Perfect (Jaccard = 1.0):    {len(perfect):>5}")
    print(f"  Good (0.8 <= J < 1.0):      {len(good):>5}")
    print(f"  Partial (0.5 <= J < 0.8):   {len(partial):>5}")
    print(f"  Poor (J < 0.5):             {len(poor):>5}")
    print(f"  Unmatched in {name1}:       {len(unmatched1):>5}")
    print(f"  Unmatched in {name2}:       {len(unmatched2):>5}")

    # Show examples of non-perfect matches
    if good:
        print(f"\nGood matches (sample):")
        for g1, m in good[:5]:
            print(f"  {g1} <-> {m['group2']}")
            print(f"    Jaccard: {m['jaccard']:.3f}, Shared: {m['shared']}, Only1: {m['only1']}, Only2: {m['only2']}")

    if partial:
        print(f"\nPartial matches (sample):")
        for g1, m in partial[:5]:
            print(f"  {g1} <-> {m['group2']}")
            print(f"    Jaccard: {m['jaccard']:.3f}, Shared: {m['shared']}, Only1: {m['only1']}, Only2: {m['only2']}")

    if poor:
        print(f"\nPoor matches (sample):")
        for g1, m in poor[:5]:
            print(f"  {g1} <-> {m['group2']}")
            print(f"    Jaccard: {m['jaccard']:.3f}, Shared: {m['shared']}, Only1: {m['only1']}, Only2: {m['only2']}")

    if unmatched1:
        print(f"\nUnmatched in {name1} (sample):")
        for g in unmatched1[:10]:
            ann = groups1[g]['annotation'][:50] if groups1[g]['annotation'] else 'N/A'
            print(f"  {g}: {len(groups1[g]['members'])} members - {ann}")

    if unmatched2:
        print(f"\nUnmatched in {name2} (sample):")
        for g in unmatched2[:10]:
            ann = groups2[g]['annotation'][:50] if groups2[g]['annotation'] else 'N/A'
            print(f"  {g}: {len(groups2[g]['members'])} members - {ann}")

    return matches, unmatched1, unmatched2


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dir1', type=Path, help='First Panaroo output directory')
    parser.add_argument('dir2', type=Path, help='Second Panaroo output directory')
    parser.add_argument('--names', nargs=2, default=None, help='Names for the two runs')

    args = parser.parse_args()

    name1 = args.names[0] if args.names else args.dir1.name
    name2 = args.names[1] if args.names else args.dir2.name

    print("=" * 60)
    print("PANAROO RUN COMPARISON (by member overlap)")
    print("=" * 60)
    print(f"\n{name1}: {args.dir1}")
    print(f"{name2}: {args.dir2}")

    # Load data
    stats1 = load_summary_stats(args.dir1)
    stats2 = load_summary_stats(args.dir2)

    groups1 = load_groups_to_members(args.dir1)
    groups2 = load_groups_to_members(args.dir2)

    # Compare
    compare_stats(stats1, stats2, name1, name2)
    match_groups_by_members(groups1, groups2, name1, name2)

    print("\n" + "=" * 60)
    print("COMPARISON COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
