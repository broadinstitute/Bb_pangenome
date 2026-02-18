#!/usr/bin/env python3
"""
Generate ENA webin-cli manifest files for genome assembly submission.

Reads:
  1. A metadata TSV (from NCBI genome submission) with columns including:
     genome_acc, biosample_accession, assembly_name, assembly_methods,
     genome_coverage, sequencing_technologies, filename
  2. An accession mapping file (TSV/CSV) that links assembly_name to
     bioproject, biosample, and SRA accessions.
  3. Assembly FASTA files.

Produces:
  - One manifest file per assembly in the output directory.
  - A summary CSV of all manifests generated.

Usage:
    python 01_generate_manifests.py \
        --metadata genome-info.tsv \
        --accessions accession_map.tsv \
        --fasta-dir /path/to/assemblies/ \
        --output-dir manifests/ \
        --study PRJNA_XXXXXXX

Notes:
    - FASTA files will be gzipped if not already compressed.
    - Assembly names in the manifest must be unique across ENA.
    - The script joins metadata and accessions on assembly_name.
    - Adjust column names in the CONFIG section if your headers differ.
"""

import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path
import json
# ============================================================
# CONFIG: adjust to match columns in metadata tables
# ============================================================

# Columns in the genome-info metadata TSV
META_ASSEMBLY_NAME = "assembly_name"
META_BIOSAMPLE = "biosample_accession"
META_COVERAGE = "genome_coverage"
META_ASSEMBLER = "assembly_methods"
META_PLATFORM = "sequencing_technologies"
META_FILENAME = "filename"

# Columns in the accession mapping file
ACC_ASSEMBLY_NAME = "isolate_name"  # join key
ACC_BIOPROJECT = "BioProject"        # PRJNA*/PRJEB*
ACC_BIOSAMPLE = "BioSample"          # SAMN*/SAMEA*
ACC_SRA = "SRA"           # SRR*/ERR*

# ENA manifest defaults
ASSEMBLY_TYPE = "isolate"
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


def ensure_gzipped(fasta_path: Path, output_dir: Path) -> Path:
    """Gzip a FASTA if not already compressed. Returns path to .gz file."""
    if fasta_path.suffix == ".gz":
        dest = output_dir / fasta_path.name
        if dest != fasta_path:
            shutil.copy2(fasta_path, dest)
        return dest

    gz_name = fasta_path.name + ".gz"
    gz_path = output_dir / gz_name
    with open(fasta_path, "rb") as f_in:
        with gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return gz_path


def find_fasta(fasta_dir: Path, filename_hint: str) -> Path | None:
    """
    Locate a FASTA file given a filename hint from the metadata.
    Searches recursively. Tries exact match first, then stem match.
    """
    # Try exact match
    exact = fasta_dir / filename_hint
    if exact.exists():
        return exact

    # Try recursive search
    hint_stem = Path(filename_hint).stem.replace(".fasta", "").replace(".fa", "").replace(".fna", "")
    for ext in (".fasta", ".fa", ".fna", ".fsa", ".fasta.gz", ".fa.gz", ".fna.gz", ".fsa.gz"): 
        matches = list(fasta_dir.rglob(f"*{hint_stem}*{ext}"))
        if len(matches) == 1:
            return matches[0]
        # If multiple matches, try tighter match
        if matches:
            tight = [m for m in matches if m.stem.replace(".fasta", "").replace(".fa", "") == hint_stem]
            if len(tight) == 1:
                return tight[0]
    return None


def write_manifest(outpath: Path, fields: dict):
    """Write a webin-cli manifest file."""
    with open(outpath, "w") as f:
        for key, value in fields.items():
            if value:
                f.write(f"{key}\t{value}\n")


def main():
    ap = argparse.ArgumentParser(description="Generate ENA webin-cli manifests.")
    ap.add_argument("--metadata", type=Path, required=True,
                    help="Genome info metadata TSV from NCBI submission")
    ap.add_argument("--accessions", type=Path, default=None,
                    help="Accession mapping file (TSV/CSV). If not provided, "
                         "biosample is taken from metadata.")
    ap.add_argument("--fasta-dir", type=Path, required=True,
                    help="Directory containing assembly FASTA files")
    ap.add_argument("--output-dir", type=Path, default=Path("manifests"),
                    help="Output directory for manifest files")
    ap.add_argument("--fasta-staging", type=Path, default=None,
                    help="Directory to stage gzipped FASTAs (default: output-dir/fastas)")
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
    fasta_staging = args.fasta_staging or (args.output_dir / "fastas")
    fasta_staging.mkdir(parents=True, exist_ok=True)

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
        filename_hint = row.get(META_FILENAME, "").strip()

        # Validate required fields
        if not biosample:
            errors.append(f"  SKIP {asm_name}: no biosample accession")
            continue

        # Find the FASTA file
        fasta_path = find_fasta(args.fasta_dir, filename_hint)
        if fasta_path is None:
            # Try assembly name as fallback
            fasta_path = find_fasta(args.fasta_dir, asm_name)
        if fasta_path is None:
            errors.append(f"  SKIP {asm_name}: FASTA not found "
                         f"(hint={filename_hint})")
            continue

        # Build unique assembly name for ENA
        ena_asm_name = f"{args.prefix}{asm_name}" if args.prefix else asm_name

        if args.dry_run:
            print(f"  [DRY RUN] {asm_name} -> {ena_asm_name}")
            print(f"    FASTA: {fasta_path}")
            print(f"    SAMPLE: {biosample} | STUDY: {bioproject}")
            continue

        gz_path = ensure_gzipped(fasta_path, fasta_staging)

        manifest_path = args.output_dir / f"manifest_{asm_name}.txt"
        manifest_fields = {
            "STUDY": bioproject,
            "SAMPLE": biosample,
            "ASSEMBLYNAME": ena_asm_name,
            "ASSEMBLY_TYPE": args.assembly_type,
            "COVERAGE": coverage,
            "PROGRAM": assembler,
            "PLATFORM": platform,
            "FASTA": str(gz_path.resolve()),
        }
        write_manifest(manifest_path, manifest_fields)

        summary.append({
            "assembly_name": asm_name,
            "ena_assembly_name": ena_asm_name,
            "biosample": biosample,
            "study": bioproject,
            "fasta": str(gz_path),
            "manifest": str(manifest_path),
            "coverage": coverage,
        })
        print(f"  [OK] {asm_name} -> {manifest_path.name}")

    print(f"\n{'='*60}")
    print(f"Generated {len(summary)} manifests in {args.output_dir}")
    if errors:
        print(f"\n{len(errors)} errors:")
        for e in errors:
            print(e)

    if summary:
        summary_path = args.output_dir / "manifest_summary.csv"
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary[0].keys())
            writer.writeheader()
            writer.writerows(summary)
        print(f"\nSummary written to {summary_path}")


if __name__ == "__main__":
    main()
