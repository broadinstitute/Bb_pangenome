#!/usr/bin/env python3
"""
Generate ENA webin-cli manifest files for annotated chromosome-level
genome assembly submission.

Reads:
  1. A metadata TSV (from NCBI genome submission) with columns including:
     genome_acc, biosample_accession, assembly_name, assembly_methods,
     genome_coverage, sequencing_technologies, filename
  2. An accession mapping file (TSV/CSV) that links assembly_name to
     bioproject, biosample, and SRA accessions.
  3. Patched EMBL flatfiles (gzipped).
  4. Chromosome list files (gzipped).
  5. Unlocalised list files (gzipped, optional per assembly).

Produces:
  - One manifest file per assembly in the output directory.
  - A summary CSV of all manifests generated.

Usage:
    python 01_generate_manifests_annotated.py \
        --metadata genome-info.tsv \
        --accessions accession_map.tsv \
        --flatfile-dir patched_embl/ \
        --chromosome-list-dir chromosome_lists/ \
        --output-dir manifests_annotated/ \
        --study PRJNA_XXXXXXX

Notes:
    - EMBL flatfiles and chromosome/unlocalised lists must be pre-gzipped.
    - Assembly names have H/P suffix stripped for ENA uniqueness.
    - The script joins metadata and accessions on assembly_name.
    - Adjust column names in the CONFIG section if headers differ.
"""

import argparse
import csv
import re
import sys
from pathlib import Path
import json

# ============================================================
# CONFIG: col names from the metadata tables
# ============================================================

# Columns in the genome-info metadata TSV
META_ASSEMBLY_NAME = "assembly_name"
META_BIOSAMPLE = "biosample_accession"
META_COVERAGE = "genome_coverage"
META_ASSEMBLER = "assembly_methods"
META_PLATFORM = "sequencing_technologies"

# Columns in the accession mapping file
ACC_ASSEMBLY_NAME = "isolate_name"  # join key
ACC_BIOPROJECT = "BioProject"        # PRJNA*/PRJEB*
ACC_BIOSAMPLE = "BioSample"          # SAMN*/SAMEA*
ACC_SRA = "SRA"                      # SRR*/ERR*

# ENA manifest defaults
ASSEMBLY_TYPE = "isolate"
MOLECULE_TYPE = "genomic DNA"
# ============================================================


def clean_field(val):
    """Convert NCBI [["value", "version"]] format to plain string."""
    try:
        parsed = json.loads(val)
        if isinstance(parsed, list):
            return " ".join(str(item) for sublist in parsed for item in sublist)
    except (json.JSONDecodeError, TypeError):
        pass
    return val


def strip_hp_suffix(name: str) -> str:
    """Strip trailing H or P from assembly name. URI88H -> URI88."""
    return re.sub(r'[HP]$', '', name)


def read_table(path: Path, delimiter=None) -> list[dict]:
    """Read a TSV or CSV, auto-detecting delimiter if not specified."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        if delimiter is None:
            if path.suffix.lower() == ".csv":
                delimiter = ","
            elif path.suffix.lower() == ".tsv":
                delimiter = "\t"
            else:
                delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader)


def write_manifest(outpath: Path, fields: dict):
    """Write a webin-cli manifest file."""
    with open(outpath, "w") as f:
        for key, value in fields.items():
            if value:
                f.write(f"{key}\t{value}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Generate ENA webin-cli manifests for annotated assembly submission."
    )
    ap.add_argument("--metadata", type=Path, required=True,
                    help="Genome info metadata TSV from NCBI submission")
    ap.add_argument("--accessions", type=Path, default=None,
                    help="Accession mapping file (TSV/CSV). If not provided, "
                         "biosample is taken from metadata.")
    ap.add_argument("--flatfile-dir", type=Path, required=True,
                    help="Directory containing gzipped patched EMBL files")
    ap.add_argument("--chromosome-list-dir", type=Path, required=True,
                    help="Directory containing gzipped chromosome/unlocalised list files")
    ap.add_argument("--output-dir", type=Path, default=Path("manifests_annotated"),
                    help="Output directory for manifest files")
    ap.add_argument("--study", type=str, required=False,
                    help="BioProject/Study accession (e.g., PRJNA_XXXXXXX)")
    ap.add_argument("--assembly-type", type=str, default=ASSEMBLY_TYPE,
                    help=f"Assembly type (default: {ASSEMBLY_TYPE})")
    ap.add_argument("--prefix", type=str, default="",
                    help="Prefix for assembly names to ensure uniqueness")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be generated without writing files")
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load metadata
    meta_rows = read_table(args.metadata)
    print(f"Loaded {len(meta_rows)} rows from metadata: {args.metadata}")
    print(f"  Columns: {list(meta_rows[0].keys()) if meta_rows else 'EMPTY'}")

    # Load accession mapping if provided
    acc_lookup = {}
    if args.accessions:
        acc_rows = read_table(args.accessions)
        print(f"Loaded {len(acc_rows)} rows from accessions: {args.accessions}")
        print(f"  Columns: {list(acc_rows[0].keys()) if acc_rows else 'EMPTY'}")
        for row in acc_rows:
            key = row.get(ACC_ASSEMBLY_NAME, "").strip()
            if key:
                acc_lookup[key] = row

    # Process each assembly
    summary = []
    errors = []

    for row in meta_rows:
        asm_name = row.get(META_ASSEMBLY_NAME, "").strip()
        if not asm_name:
            continue

        # Get biosample - prefer accession mapping, fall back to metadata
        if asm_name in acc_lookup:
            biosample = acc_lookup[asm_name].get(ACC_BIOSAMPLE, "").strip()
            bioproject = acc_lookup[asm_name].get(ACC_BIOPROJECT, args.study).strip()
        else:
            biosample = row.get(META_BIOSAMPLE, "").strip()
            bioproject = args.study

        coverage = row.get(META_COVERAGE, "").strip()
        coverage = coverage.rstrip("xX")
        assembler = clean_field(row.get(META_ASSEMBLER, "").strip())
        platform = clean_field(row.get(META_PLATFORM, "").strip())

        # Validate required fields
        if not biosample:
            errors.append(f"  SKIP {asm_name}: no biosample accession")
            continue

        # Find EMBL flatfile
        embl_path = args.flatfile_dir / f"{asm_name}.embl.gz"
        if not embl_path.exists():
            errors.append(f"  SKIP {asm_name}: EMBL not found at {embl_path}")
            continue

        # Find chromosome list (required for chromosome-level)
        chrom_path = args.chromosome_list_dir / f"{asm_name}.chromosome_list.tsv.gz"
        if not chrom_path.exists():
            errors.append(f"  SKIP {asm_name}: chromosome list not found at {chrom_path}")
            continue

        # Find unlocalised list (optional)
        unloc_path = args.chromosome_list_dir / f"{asm_name}.unlocalised_list.tsv.gz"
        has_unloc = unloc_path.exists()

        # Build unique assembly name: strip H/P suffix
        #ena_asm_name = strip_hp_suffix(asm_name)
        ena_asm_name = f"{asm_name}_ann_v1"
        if args.prefix:
            ena_asm_name = f"{args.prefix}{ena_asm_name}"

        if args.dry_run:
            print(f"  [DRY RUN] {asm_name} -> {ena_asm_name}")
            print(f"    FLATFILE: {embl_path}")
            print(f"    CHROMOSOME_LIST: {chrom_path}")
            if has_unloc:
                print(f"    UNLOCALISED_LIST: {unloc_path}")
            print(f"    SAMPLE: {biosample} | STUDY: {bioproject}")
            continue

        # Write manifest
        manifest_path = args.output_dir / f"manifest_{asm_name}.txt"
        manifest_fields = {
            "STUDY": bioproject,
            "SAMPLE": biosample,
            "ASSEMBLYNAME": ena_asm_name,
            "ASSEMBLY_TYPE": args.assembly_type,
            "COVERAGE": coverage,
            "PROGRAM": assembler,
            "PLATFORM": platform,
            "MOLECULETYPE": MOLECULE_TYPE,
            "FLATFILE": str(embl_path.resolve()),
            "CHROMOSOME_LIST": str(chrom_path.resolve()),
        }
        if has_unloc:
            manifest_fields["UNLOCALISED_LIST"] = str(unloc_path.resolve())

        write_manifest(manifest_path, manifest_fields)

        summary.append({
            "assembly_name": asm_name,
            "ena_assembly_name": ena_asm_name,
            "biosample": biosample,
            "study": bioproject,
            "flatfile": str(embl_path),
            "chromosome_list": str(chrom_path),
            "unlocalised_list": str(unloc_path) if has_unloc else "",
            "manifest": str(manifest_path),
            "coverage": coverage,
        })
        print(f"  [OK] {asm_name} -> {manifest_path.name}"
              f"{' (+unlocalised)' if has_unloc else ''}")

    # Report
    print(f"\n{'='*60}")
    print(f"Generated {len(summary)} manifests in {args.output_dir}")
    if errors:
        print(f"\n{len(errors)} errors:")
        for e in errors:
            print(e)

    # Write summary
    if summary:
        summary_path = args.output_dir / "manifest_summary.csv"
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary[0].keys())
            writer.writeheader()
            writer.writerows(summary)
        print(f"\nSummary written to {summary_path}")


if __name__ == "__main__":
    main()
