#!/usr/bin/env python3
"""
Replace locus tags in EMBL annotation files using a mapping table.

Reads EMBL files from assemblies/{isolate}/{isolate}.embl,
replaces old (bakta) locus tags with new (NCBI-assigned) tags,
and writes updated files to the output directory.

Usage:
    python replace_locus_tags.py \
        --assemblies /path/to/assemblies/ \
        --map locus_tag_map.csv \
        --output /path/to/output/

Map file format (CSV/TSV with header):
    isolate,old,new
    B418P,BAKTA_001,LOC_001
    B331P,BAKTA_002,LOC_002
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path


def read_map(map_path: Path) -> list[dict]:
    """Read locus tag map, auto-detecting delimiter."""
    with open(map_path, newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        if map_path.suffix.lower() == ".tsv":
            delimiter = "\t"
        elif map_path.suffix.lower() == ".csv":
            delimiter = ","
        else:
            delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader)


def replace_in_embl(embl_path: Path, old_tag: str, new_tag: str) -> tuple[str, int]:
    """Replace all occurrences of old locus tag with new in EMBL file.
    Returns the modified text and the number of replacements made."""
    text = embl_path.read_text(encoding="utf-8")
    count = text.count(old_tag)
    updated = text.replace(old_tag, new_tag)
    return updated, count


def main():
    ap = argparse.ArgumentParser(
        description="Replace locus tags in EMBL files using a mapping table."
    )
    ap.add_argument("--assemblies", type=Path, required=True,
                    help="Input assemblies directory (structure: {isolate}/{isolate}.embl)")
    ap.add_argument("--map", type=Path, required=True,
                    help="Locus tag mapping file (columns: isolate, old, new)")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output directory for updated EMBL files")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would change without writing files")
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    # Read mapping
    tag_map = read_map(args.map)
    print(f"Loaded {len(tag_map)} entries from {args.map}")

    success = 0
    skipped = 0
    errors = []

    for row in tag_map:
        isolate = row.get("isolate", "").strip()
        old_tag = row.get("old", "").strip()
        new_tag = row.get("new", "").strip()

        if not all([isolate, old_tag, new_tag]):
            errors.append(f"  SKIP: incomplete row: {row}")
            skipped += 1
            continue

        # Find EMBL file
        embl_path = args.assemblies / isolate / f"{isolate}.embl"
        if not embl_path.exists():
            # Try without subdirectory
            embl_path = args.assemblies / f"{isolate}.embl"
        if not embl_path.exists():
            errors.append(f"  SKIP {isolate}: EMBL file not found")
            skipped += 1
            continue

        updated_text, count = replace_in_embl(embl_path, old_tag, new_tag)

        if count == 0:
            errors.append(f"  WARN {isolate}: no occurrences of '{old_tag}' found")

        if args.dry_run:
            print(f"  [DRY RUN] {isolate}: {old_tag} -> {new_tag} ({count} replacements)")
            continue

        out_path = args.output / f"{isolate}.embl"
        out_path.write_text(updated_text, encoding="utf-8")
        print(f"  [OK] {isolate}: {old_tag} -> {new_tag} ({count} replacements) -> {out_path.name}")
        success += 1

    print(f"\n{'='*60}")
    print(f"Done: {success} updated, {skipped} skipped")
    if errors:
        print(f"\nIssues:")
        for e in errors:
            print(e)


if __name__ == "__main__":
    main()
