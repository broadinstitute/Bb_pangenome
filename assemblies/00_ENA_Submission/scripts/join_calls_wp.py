#!/usr/bin/env python3
"""
Join existing replicon calls (10.2) with whole-plasmid (wp) alignment stats
for chromosome list generation.

Takes the replicon name from the calls file and the alignment stats
(ref_length, ref_covered_length, query_coverage_percent, overall_percent_identity)
from the wp best hits table.

Usage:
    python join_calls_wp.py \
        --calls best_hits_1000bp_v10.2.csv \
        --wp calls/wp/tables/summary_best_hits.tsv \
        --output calls_with_wp_stats.tsv
"""

import argparse
import pandas as pd
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(
        description="Join replicon calls with wp alignment stats."
    )
    ap.add_argument("--calls", type=Path, required=True,
                    help="Published calls file (e.g., best_hits_1000bp_v10.2.csv)")
    ap.add_argument("--wp", type=Path, required=True,
                    help="WP best hits table with alignment stats")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output merged table")
    ap.add_argument("--calls-contig-col", default="contig_id",
                    help="Contig ID column in calls file (default: contig_id)")
    ap.add_argument("--calls-assembly-col", default="assembly_id",
                    help="Assembly ID column in calls file (default: assembly_id)")
    ap.add_argument("--calls-call-col", default="plasmid_name",
                    help="Call column in calls file (default: plasmid_name)")
    ap.add_argument("--wp-contig-col", default="contig_id",
                    help="Contig ID column in wp file (default: contig_id)")
    ap.add_argument("--wp-assembly-col", default="assembly_id",
                    help="Assembly ID column in wp file (default: assembly_id)")
    args = ap.parse_args()

    calls_sep = "," if args.calls.suffix == ".csv" else "\t"
    calls = pd.read_csv(args.calls, sep=calls_sep, encoding="utf-8-sig")
    print(f"Calls: {len(calls)} rows from {args.calls.name}")
    wp_sep = "," if args.wp.suffix == ".csv" else "\t"
    wp = pd.read_csv(args.wp, sep=wp_sep, encoding="utf-8-sig")
    print(f"WP hits: {len(wp)} rows from {args.wp.name}")

    import re
    def strip_brackets(s):
        return re.sub(r'\s*\[.*?\]', '', str(s)).strip()

    calls["_join_contig"] = calls[args.calls_contig_col].apply(strip_brackets)
    calls["_join_assembly"] = calls[args.calls_assembly_col].str.strip()
    wp["_join_contig"] = wp[args.wp_contig_col].apply(strip_brackets)
    wp["_join_assembly"] = wp[args.wp_assembly_col].str.strip()
    wp_stats = wp[["_join_assembly", "_join_contig",
                   "ref_length_wp", "ref_covered_length_wp",
                   "query_coverage_percent_wp", "overall_percent_identity_wp",
                   "plasmid_name_wp"]].copy()
    wp_stats = wp_stats.rename(columns={
        "ref_length_wp": "wp_ref_length",
        "ref_covered_length_wp": "wp_ref_covered_length",
        "query_coverage_percent_wp": "wp_query_coverage_percent",
        "overall_percent_identity_wp": "wp_overall_percent_identity",
        "plasmid_name_wp": "wp_plasmid_name",
    })
    merged = calls.merge(
        wp_stats,
        on=["_join_assembly", "_join_contig"],
        how="left"
    )
    has_wp = merged["wp_ref_length"].notna().sum()
    no_wp = merged["wp_ref_length"].isna().sum()
    print(f"\nJoin results:")
    print(f"  With wp stats: {has_wp}")
    print(f"  Without wp stats: {no_wp}")

    out = merged[[
        args.calls_assembly_col, args.calls_contig_col, "contig_len",
        args.calls_call_col,
        "wp_plasmid_name",
        "wp_ref_length", "wp_ref_covered_length",
        "wp_query_coverage_percent", "wp_overall_percent_identity",
    ]].copy()

    out = out.rename(columns={
        args.calls_call_col: "plasmid_name",
        "wp_ref_length": "ref_length",
        "wp_ref_covered_length": "ref_covered_length",
        "wp_query_coverage_percent": "query_coverage_percent",
        "wp_overall_percent_identity": "overall_percent_identity",
    })

    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {len(out)} rows -> {args.output}")


if __name__ == "__main__":
    main()
