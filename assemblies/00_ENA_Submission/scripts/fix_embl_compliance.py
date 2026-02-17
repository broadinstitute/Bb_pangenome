#!/usr/bin/env python3
"""
Fix EMBL flatfile compliance issues for ENA submission.

Fixes:
  1. Removes /tag_peptide qualifiers (not accepted by ENA validator)
  2. Corrects ID line topology (linear/circular) to match chromosome list

Reads chromosome list files to determine which contigs should be circular.
Processes all EMBL files in the input directory.

Usage:
    python fix_embl_compliance.py \
        --input-dir patched_embl/ \
        --output-dir fixed_embl/ \
        --chromosome-list-dir chromosome_lists/

    # In-place fix (overwrites input files):
    python fix_embl_compliance.py \
        --input-dir patched_embl/ \
        --output-dir patched_embl/ \
        --chromosome-list-dir chromosome_lists/
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def load_circular_contigs(chromosome_list_dir: Path) -> dict[str, set[str]]:
    """Load chromosome lists and return {assembly_id: set of circular contig names}.
    
    Reads both .tsv and .tsv.gz files. Identifies circular contigs from
    the CHROMOSOME_TYPE column (third column containing 'Circular').
    """
    circular = defaultdict(set)

    for chrom_file in sorted(chromosome_list_dir.glob("*.chromosome_list.tsv")):
        assembly_id = chrom_file.name.replace(".chromosome_list.tsv", "")
        with open(chrom_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    object_name = parts[0].strip()
                    chrom_type = parts[2].strip()
                    if "Circular" in chrom_type:
                        circular[assembly_id].add(object_name)

    return dict(circular)


def fix_embl_file(input_path: Path, output_path: Path,
                  circular_contigs: set[str],
                  dry_run: bool = False) -> dict:
    """Fix a single EMBL file.
    
    Returns dict with counts of fixes applied.
    """
    with open(input_path, "r") as f:
        lines = f.readlines()

    fixed = []
    stats = {
        "tag_peptide_removed": 0,
        "topology_fixed": 0,
        "entries": 0,
    } # tag_peptide needs to be dropped from our data explicitly.

    current_contig = None
    i = 0
    while i < len(lines):
        line = lines[i]

        # Track current entry by ID line
        # ID   contig_1; SV 1; linear; genomic DNA; STD; UNK; 903099 BP.
        if line.startswith("ID   "):
            stats["entries"] += 1
            # Extract contig name from ID line
            id_match = re.match(r'ID\s+(\S+?)\s*;', line)
            if id_match:
                current_contig = id_match.group(1)

            # Fix topology if this contig should be circular
            if current_contig and current_contig in circular_contigs:
                if "linear" in line.lower() and "circular" not in line.lower():
                    line = line.replace("linear", "circular").replace("LINEAR", "CIRCULAR")
                    stats["topology_fixed"] += 1

            fixed.append(line)
            i += 1
            continue

        # Remove /tag_peptide qualifier lines
        # FT                   /tag_peptide="complement(46644..46691)"
        if line.startswith("FT") and "/tag_peptide=" in line:
            stats["tag_peptide_removed"] += 1
            i += 1
            # Skip continuation lines (FT lines starting with spaces after qualifier)
            while i < len(lines) and lines[i].startswith("FT") and \
                  re.match(r'^FT\s{19,}[^/]', lines[i]) and \
                  "/tag_peptide" not in lines[i - 1]:
                # This handles multi-line qualifier values
                i += 1
            continue

        # Also remove the tmRNA feature's tag_peptide — check if it's a
        # standalone FT qualifier line
        fixed.append(line)
        i += 1

    if not dry_run:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.writelines(fixed)

    return stats


def main():
    ap = argparse.ArgumentParser(
        description="Fix EMBL flatfile compliance issues for ENA submission."
    )
    ap.add_argument("--input-dir", type=Path, required=True,
                    help="Directory containing .embl files to fix")
    ap.add_argument("--output-dir", type=Path, required=True,
                    help="Output directory for fixed .embl files")
    ap.add_argument("--chromosome-list-dir", type=Path, required=True,
                    help="Directory containing chromosome list TSV files")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be changed without writing files")
    args = ap.parse_args()

    # Load circular contig sets from chromosome lists
    circular_map = load_circular_contigs(args.chromosome_list_dir)
    total_circular = sum(len(v) for v in circular_map.values())
    print(f"Loaded chromosome lists for {len(circular_map)} assemblies")
    print(f"Total circular contigs identified: {total_circular}")

    # Process EMBL files
    embl_files = sorted(args.input_dir.glob("*.embl"))
    print(f"Found {len(embl_files)} EMBL files in {args.input_dir}\n")

    if not embl_files:
        print("No .embl files found!")
        sys.exit(1)

    total_stats = {
        "files": 0,
        "entries": 0,
        "tag_peptide_removed": 0,
        "topology_fixed": 0,
    }

    for embl_path in embl_files:
        assembly_id = embl_path.stem
        circular_contigs = circular_map.get(assembly_id, set())
        output_path = args.output_dir / embl_path.name

        stats = fix_embl_file(
            embl_path, output_path,
            circular_contigs=circular_contigs,
            dry_run=args.dry_run,
        )

        prefix = "[DRY RUN] " if args.dry_run else "[OK] "
        fixes = []
        if stats["tag_peptide_removed"]:
            fixes.append(f"{stats['tag_peptide_removed']} tag_peptide")
        if stats["topology_fixed"]:
            fixes.append(f"{stats['topology_fixed']} topology")

        fix_str = ", ".join(fixes) if fixes else "no fixes needed"
        print(f"  {prefix}{embl_path.name}: {stats['entries']} entries — {fix_str}")

        total_stats["files"] += 1
        for k in ["entries", "tag_peptide_removed", "topology_fixed"]:
            total_stats[k] += stats[k]

    print(f"\n{'='*60}")
    print(f"Files processed: {total_stats['files']}")
    print(f"Total entries: {total_stats['entries']}")
    print(f"tag_peptide qualifiers removed: {total_stats['tag_peptide_removed']}")
    print(f"Topology conflicts fixed: {total_stats['topology_fixed']}")

    if args.dry_run:
        print("\n(dry run — no files written)")


if __name__ == "__main__":
    main()
