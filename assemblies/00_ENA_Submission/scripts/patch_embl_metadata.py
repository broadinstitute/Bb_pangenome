#!/usr/bin/env python3
"""
Patch EMBL flatfiles with organism metadata for ENA submission.

For each entry in a multi-entry EMBL file, inserts:
  1. OS line: organism species name
  2. OC lines: full taxonomy lineage
  3. FT source qualifiers: /organism and /strain

Usage:
    python patch_embl_metadata.py \
        --input-dir cleaned_embl/ \
        --output-dir patched_embl/ \
        --strain-map strain_map.tsv \
        --organism "Borrelia burgdorferi" \
        --taxonomy "Bacteria; Spirochaetota; Spirochaetia; Spirochaetales; Borreliaceae; Borrelia;"

    strain_map.tsv format (tab-separated, no header):
        assembly_id<TAB>strain_name
        URI88H<TAB>URI88H
        B500P<TAB>B500P
"""

import argparse
import re
import sys
from pathlib import Path


def load_strain_map(path: Path) -> dict:
    """Load assembly_id -> strain_name mapping."""
    mapping = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mapping[parts[0].strip()] = parts[1].strip()
            elif len(parts) == 1:
                # assembly_id is the strain name
                mapping[parts[0].strip()] = parts[0].strip()
    return mapping


def patch_embl_file(input_path: Path, output_path: Path,
                    organism: str, taxonomy: str, strain: str,
                    dry_run: bool = False) -> dict:
    """Patch a single EMBL file with metadata.
    
    Returns dict with counts of patches applied.
    """
    with open(input_path, "r") as f:
        lines = f.readlines()

    patched = []
    stats = {"entries": 0, "os_patched": 0, "oc_patched": 0,
             "organism_added": 0, "strain_added": 0}

    # Format taxonomy into OC lines (max 80 chars per line)
    oc_lines = []
    tax_parts = [t.strip() for t in taxonomy.rstrip(";").split(";")]
    current_line = "OC   "
    for i, part in enumerate(tax_parts):
        suffix = "; " if i < len(tax_parts) - 1 else ";"
        candidate = current_line + part + suffix
        if len(candidate) > 78 and current_line != "OC   ":
            oc_lines.append(current_line.rstrip() + "\n")
            current_line = "OC   " + part + suffix
        else:
            current_line = candidate
    if current_line.strip() != "OC":
        oc_lines.append(current_line.rstrip() + "\n")

    i = 0
    while i < len(lines):
        line = lines[i]

        # Patch OS line: "OS   ." -> "OS   Borrelia burgdorferi"
        if line.startswith("OS   ") and line.strip() == "OS   .":
            patched.append(f"OS   {organism}\n")
            stats["os_patched"] += 1
            i += 1
            continue

        # Patch OC line: "OC   ." -> full taxonomy
        if line.startswith("OC   ") and line.strip() == "OC   .":
            for oc_line in oc_lines:
                patched.append(oc_line)
            stats["oc_patched"] += 1
            i += 1
            continue

        # After /mol_type line in FT source, insert /organism and /strain
        if (line.startswith("FT") and '/mol_type="genomic DNA"' in line):
            patched.append(line)
            # Check if /organism already present in next lines
            has_organism = False
            has_strain = False
            j = i + 1
            while j < len(lines) and lines[j].startswith("FT"):
                if "/organism=" in lines[j]:
                    has_organism = True
                if "/strain=" in lines[j]:
                    has_strain = True
                if lines[j].startswith("FT   ") and not lines[j].startswith("FT                   "):
                    # Next feature, stop checking
                    break
                j += 1

            if not has_organism:
                patched.append(f'FT                   /organism="{organism}"\n')
                stats["organism_added"] += 1
            if not has_strain:
                patched.append(f'FT                   /strain="{strain}"\n')
                stats["strain_added"] += 1
            i += 1
            continue

        # Count entries
        if line.startswith("//"):
            stats["entries"] += 1

        patched.append(line)
        i += 1

    if not dry_run:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.writelines(patched)

    return stats


def main():
    ap = argparse.ArgumentParser(
        description="Patch EMBL flatfiles with organism metadata for ENA submission."
    )
    ap.add_argument("--input-dir", type=Path, required=True,
                    help="Directory containing .embl files")
    ap.add_argument("--output-dir", type=Path, required=True,
                    help="Output directory for patched .embl files")
    ap.add_argument("--strain-map", type=Path, required=True,
                    help="TSV mapping assembly_id to strain name (no header)")
    ap.add_argument("--organism", type=str, required=True,
                    help='Organism name, e.g. "Borrelia burgdorferi"')
    ap.add_argument("--taxonomy", type=str, required=True,
                    help='Semicolon-separated taxonomy, e.g. '
                         '"Bacteria; Spirochaetota; Spirochaetia; ..."')
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be changed without writing files")
    args = ap.parse_args()

    strain_map = load_strain_map(args.strain_map)
    print(f"Loaded {len(strain_map)} strain mappings")

    embl_files = sorted(args.input_dir.glob("*.embl"))
    print(f"Found {len(embl_files)} EMBL files in {args.input_dir}")

    if not embl_files:
        print("No .embl files found!")
        sys.exit(1)

    total_stats = {"files": 0, "entries": 0, "os_patched": 0,
                   "oc_patched": 0, "organism_added": 0, "strain_added": 0}
    missing_strains = []

    for embl_path in embl_files:
        # Extract assembly_id from filename (e.g., URI88H.embl -> URI88H)
        assembly_id = embl_path.stem

        if assembly_id not in strain_map:
            missing_strains.append(assembly_id)
            strain = assembly_id  # fallback: use assembly_id as strain
        else:
            strain = strain_map[assembly_id]

        output_path = args.output_dir / embl_path.name

        stats = patch_embl_file(
            embl_path, output_path,
            organism=args.organism,
            taxonomy=args.taxonomy,
            strain=strain,
            dry_run=args.dry_run,
        )

        prefix = "[DRY RUN] " if args.dry_run else "[OK] "
        print(f"  {prefix}{embl_path.name}: {stats['entries']} entries, "
              f"{stats['os_patched']} OS, {stats['oc_patched']} OC, "
              f"{stats['organism_added']} /organism, {stats['strain_added']} /strain")

        total_stats["files"] += 1
        for k in ["entries", "os_patched", "oc_patched",
                   "organism_added", "strain_added"]:
            total_stats[k] += stats[k]

    print(f"\n{'='*60}")
    print(f"Files processed: {total_stats['files']}")
    print(f"Total entries: {total_stats['entries']}")
    print(f"OS lines patched: {total_stats['os_patched']}")
    print(f"OC lines patched: {total_stats['oc_patched']}")
    print(f"/organism added: {total_stats['organism_added']}")
    print(f"/strain added: {total_stats['strain_added']}")

    if missing_strains:
        print(f"\n⚠  {len(missing_strains)} assemblies not in strain map "
              f"(used assembly_id as strain):")
        for s in missing_strains:
            print(f"  {s}")


if __name__ == "__main__":
    main()
