#!/usr/bin/env python3
"""
Generate ENA chromosome list and unlocalised list files from v11 calls.

For each assembly, produces:
  1. chromosome_list.tsv  — one entry per replicon (longest contig wins)
  2. unlocalised_list.tsv — fragment contigs associated with a replicon
                            but not the primary sequence

Two modes:
  --mode complete   Only include contigs passing completeness thresholds
                    in the chromosome list (strict, default)
  --mode classified All classified contigs get placed; longest per replicon
                    in chromosome list, rest in unlocalised list

Chromosome list format (ENA):
    OBJECT_NAME  CHROMOSOME_NAME  CHROMOSOME_TYPE

Unlocalised list format (ENA):
    OBJECT_NAME  CHROMOSOME_NAME

Usage:
    python generate_chromosome_lists.py \
        --classifier v11_with_wp_stats.tsv \
        --output-dir chromosome_lists/ \
        --mode classified \
        --dry-run
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from collections import defaultdict


DEFAULT_REF_COV = 0.95
DEFAULT_QUERY_COV = 0.90
DEFAULT_IDENTITY = 0.90


def parse_contig_name(contig_id: str) -> str:
    """Extract the base contig name, stripping bracket metadata.
    e.g., 'contig_1 [gcode=11] [topology=linear]' -> 'contig_1'
    """
    return contig_id.strip().split()[0]


def infer_topology(plasmid_name: str) -> str:
    """Infer topology from replicon name."""
    if not isinstance(plasmid_name, str):
        return "linear"
    n = plasmid_name.lower().strip()
    if n.startswith("cp"):
        return "circular"
    if n.startswith("lp") or n == "chromosome":
        return "linear"
    return "linear"



def sanitize_replicon_name(name: str) -> str:
    """Make replicon name ENA-compliant.
    - 'chromosome' -> 'main' (reserved word)
    - '+' -> '-' (invalid character)
    """
    name = name.replace('+', '-')
    if name.lower() == 'chromosome':
        name = 'main'
    return name


def infer_chromosome_type(plasmid_name: str) -> str:
    """Determine if this is a Chromosome or Plasmid for ENA."""
    if not isinstance(plasmid_name, str):
        return "Plasmid"
    if plasmid_name.lower().strip() == "chromosome":
        return "Chromosome"
    return "Plasmid"


def format_topology(topology: str, chrom_type: str) -> str:
    """Format topology + type for ENA chromosome list.
    ENA expects: Linear-Chromosome, Circular-Plasmid, etc.
    """
    return f"{topology.capitalize()}-{chrom_type}"


def compute_ref_coverage(ref_covered_length, ref_length) -> float | None:
    """Compute fraction of reference covered."""
    try:
        rcl = float(ref_covered_length)
        rl = float(ref_length)
        if rl <= 0:
            return None
        return rcl / rl
    except (ValueError, TypeError):
        return None


def passes_completeness(row, ref_cov_thresh, query_cov_thresh, identity_thresh) -> bool:
    """Check if a contig passes completeness thresholds."""
    ref_cov = compute_ref_coverage(
        row.get("ref_covered_length", ""),
        row.get("ref_length", "")
    )
    try:
        q_cov = float(row.get("query_coverage_percent", 0)) / 100.0
    except (ValueError, TypeError):
        q_cov = 0.0
    try:
        ident = float(row.get("overall_percent_identity", 0)) / 100.0
    except (ValueError, TypeError):
        ident = 0.0

    return (ref_cov is not None
            and ref_cov >= ref_cov_thresh
            and q_cov >= query_cov_thresh
            and ident >= identity_thresh)


def main():
    ap = argparse.ArgumentParser(
        description="Generate ENA chromosome list and unlocalised list files."
    )
    ap.add_argument("--classifier", type=Path, required=True,
                    help="Classifier output with alignment stats")
    ap.add_argument("--output-dir", type=Path, required=True,
                    help="Output directory for chromosome/unlocalised list files")
    ap.add_argument("--mode", choices=["complete", "classified"], default="complete",
                    help="'complete' = strict completeness thresholds; "
                         "'classified' = all classified contigs, longest per replicon "
                         "(default: complete)")
    ap.add_argument("--ref-cov", type=float, default=DEFAULT_REF_COV,
                    help=f"Min reference coverage (complete mode only, default: {DEFAULT_REF_COV})")
    ap.add_argument("--query-cov", type=float, default=DEFAULT_QUERY_COV,
                    help=f"Min query coverage (complete mode only, default: {DEFAULT_QUERY_COV})")
    ap.add_argument("--identity", type=float, default=DEFAULT_IDENTITY,
                    help=f"Min identity (complete mode only, default: {DEFAULT_IDENTITY})")
    ap.add_argument("--call-col", default="plasmid_name",
                    help="Column name for replicon call (default: plasmid_name)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be generated without writing files")
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Read classifier output
    with open(args.classifier, newline="", encoding="utf-8-sig") as f:
        delimiter = "," if args.classifier.suffix.lower() == ".csv" else "\t"
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = list(reader)

    print(f"Loaded {len(rows)} rows from {args.classifier}")
    print(f"Mode: {args.mode}")
    if args.mode == "complete":
        print(f"Thresholds: ref_cov={args.ref_cov:.0%}, "
              f"query_cov={args.query_cov:.0%}, identity={args.identity:.0%}")

    # Group by assembly
    assemblies = defaultdict(list)
    for row in rows:
        aid = row.get("assembly_id", "").strip()
        if aid:
            assemblies[aid].append(row)

    print(f"Found {len(assemblies)} assemblies\n")

    total_chrom_files = 0
    total_unloc_files = 0
    total_placed = 0
    total_unlocalised = 0
    total_unplaced = 0
    assemblies_without_entries = []

    for assembly_id in sorted(assemblies.keys()):
        contigs = assemblies[assembly_id]
        chrom_lines = []
        unloc_lines = []

        if args.mode == "classified":
            # Group contigs by replicon call
            by_replicon = defaultdict(list)
            for row in contigs:
                call = row.get(args.call_col, "").strip()
                if not call or call.lower() == "unclassified":
                    total_unplaced += 1
                    continue
                try:
                    clen = int(float(row.get("contig_len", 0)))
                except (ValueError, TypeError):
                    clen = 0
                by_replicon[call].append((clen, row))

            for replicon, contig_list in sorted(by_replicon.items()):
                # Sort by length descending — longest is the primary
                contig_list.sort(key=lambda x: x[0], reverse=True)

                ena_name = sanitize_replicon_name(replicon)

                for i, (clen, row) in enumerate(contig_list):
                    contig_id = row.get("contig_id", "").strip()
                    object_name = parse_contig_name(contig_id)
                    topology = infer_topology(replicon)
                    chrom_type = infer_chromosome_type(replicon)
                    topo_type = format_topology(topology, chrom_type)

                    if i == 0:
                        # Primary — chromosome list
                        chrom_lines.append(
                            f"{object_name}\t{ena_name}\t{topo_type}"
                        )
                        total_placed += 1
                    else:
                        # Fragment — unlocalised list
                        unloc_lines.append(
                            f"{object_name}\t{ena_name}"
                        )
                        total_unlocalised += 1

        elif args.mode == "complete":
            for row in contigs:
                contig_id = row.get("contig_id", "").strip()
                plasmid_name = row.get(args.call_col, "").strip()

                if not contig_id or not plasmid_name:
                    continue

                if passes_completeness(row, args.ref_cov, args.query_cov, args.identity):
                    object_name = parse_contig_name(contig_id)
                    topology = infer_topology(plasmid_name)
                    chrom_type = infer_chromosome_type(plasmid_name)
                    topo_type = format_topology(topology, chrom_type)
                    ena_name = sanitize_replicon_name(plasmid_name)
                    chrom_lines.append(
                        f"{object_name}\t{ena_name}\t{topo_type}"
                    )
                    total_placed += 1
                else:
                    total_unplaced += 1

        if not chrom_lines:
            assemblies_without_entries.append(assembly_id)
            continue

        if args.dry_run:
            print(f"  [DRY RUN] {assembly_id}: "
                  f"{len(chrom_lines)} placed, {len(unloc_lines)} unlocalised")
            for line in chrom_lines:
                print(f"    [CHROM] {line}")
            for line in unloc_lines:
                print(f"    [UNLOC] {line}")
            continue

        # Write chromosome list
        chrom_path = args.output_dir / f"{assembly_id}.chromosome_list.tsv"
        with open(chrom_path, "w") as f:
            for line in chrom_lines:
                f.write(line + "\n")
        total_chrom_files += 1

        # Write unlocalised list (only if there are unlocalised contigs)
        if unloc_lines:
            unloc_path = args.output_dir / f"{assembly_id}.unlocalised_list.tsv"
            with open(unloc_path, "w") as f:
                for line in unloc_lines:
                    f.write(line + "\n")
            total_unloc_files += 1

        print(f"  [OK] {assembly_id}: {len(chrom_lines)} placed, "
              f"{len(unloc_lines)} unlocalised")

    # Summary
    print(f"\n{'='*60}")
    print(f"Generated {total_chrom_files} chromosome list files")
    print(f"Generated {total_unloc_files} unlocalised list files")
    print(f"Placed replicons (chromosome list): {total_placed}")
    print(f"Unlocalised fragments: {total_unlocalised}")
    print(f"Unplaced/unclassified: {total_unplaced}")
    if assemblies_without_entries:
        print(f"\n{len(assemblies_without_entries)} assemblies with NO entries "
              f"(will remain contig-level):")
        for a in assemblies_without_entries:
            print(f"  {a}")


if __name__ == "__main__":
    main()
