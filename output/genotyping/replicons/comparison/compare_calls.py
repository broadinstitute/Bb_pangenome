#!/usr/bin/env python3
"""
Compare existing and new classifier calls, categorize differences,
and produce a homogenized call table.

Usage:
    python compare_calls.py \
        --old old_summary.tsv \
        --new new_summary.tsv \
        --old-contig-col contig_id \
        --old-call-col final_call \
        --new-contig-col contig_id \
        --new-call-col final_call \
        --old-assembly-col assembly_id \
        --new-assembly-col assembly_id \
        --output comparison.tsv

Categories:
    exact_match        - identical calls
    extra_annotation   - new call contains old call but with additions (e.g. lp28-4 -> lp28-4:::lp17) # This is to be dealt with in detail in future work.
    base_match         - old call found within new call after stripping annotations
    different          - genuinely different classification
    old_only           - contig in old but not new
    new_only           - contig in new but not old
    new_unclassified   - new call is NA/empty but old had a call 
    old_unclassified   - old call is NA/empty but new has a call # due to length filtering of old.
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from collections import Counter


def normalize_call(call: str) -> str:
    """Normalize a call string for comparison."""
    if not call or call.lower() in ("na", "nan", "none", "", "unclassified"):
        return ""
    return call.strip().lower()


def is_empty(call: str) -> bool:
    return normalize_call(call) == ""


def extract_base_plasmid_calls(call: str) -> set:
    """Extract individual replicon names from a potential compound call.
    e.g. 'lp28-4:::lp17' -> {'lp28-4', 'lp17'}
    e.g. 'lp28-4+lp28-3' -> {'lp28-4+lp28-3'} (fusion, keep as unit, this comes from the DB.)
    """
    if is_empty(call):
        return set()
    # Split on ::: (multi-locus separator) but NOT on + (fusion)
    # + is from the db and ::: is a feature I introduced into the plasmid caller 
    # after initial submissiond
    parts = re.split(r':::', normalize_call(call))
    return {p.strip() for p in parts if p.strip()}


def strip_suffix(call: str) -> str:
    """Strip trailing annotation markers like * from a call.
    e.g. 'lp28-1*' -> 'lp28-1'
    """
    return re.sub(r'\*+$', '', call.strip())


def get_replicon_family(call: str) -> str | None:
    """Extract the replicon family prefix.
    e.g. 'cp32-12' -> 'cp32', 'lp28-4' -> 'lp28', 'chromosome' -> 'chromosome'
    """
    call = normalize_call(strip_suffix(call))
    if not call or call == "unclassified":
        return None
    # cp32-12 -> cp32, lp28-4 -> lp28, lp54 -> lp54, chromosome -> chromosome
    m = re.match(r'^((?:cp|lp)\d+)', call)
    if m:
        return m.group(1)
    return call


def categorize(old_call: str, new_call: str) -> tuple[str, str]:
    """
    Categorize the relationship between old and new calls.
    Returns (category, resolution).
    Resolution is the recommended homogenized call.
    """
    old_norm = normalize_call(old_call)
    new_norm = normalize_call(new_call)

    # Both empty
    if is_empty(old_call) and is_empty(new_call):
        return "exact_match", ""

    # New lost classification
    if not is_empty(old_call) and is_empty(new_call):
        return "new_unclassified", old_call.strip()

    # Old was unclassified, new has call
    if is_empty(old_call) and not is_empty(new_call):
        return "old_unclassified", new_call.strip()

    # Exact match
    if old_norm == new_norm:
        return "exact_match", old_call.strip()

    # Annotation suffix: lp28-1 vs lp28-1* (same base, trailing *)
    if strip_suffix(old_norm) == strip_suffix(new_norm):
        return "annotation_suffix", old_call.strip()

    # Check if old call is contained within new compound call
    old_base_calls = extract_base_plasmid_calls(old_call)
    new_base_calls = extract_base_plasmid_calls(new_call)

    if old_base_plasmid_calls and old_base_plasmid_calls.issubset(new_base_plasmid_calls):
        # Old call is a subset of new — new has extra annotations
        return "extra_annotation", old_call.strip()

    if new_base_plasmid_calls and new_base_plasmid_calls.issubset(old_base_plasmid_calls):
        # New is subset of old — old had more info
        return "base_match", old_call.strip()

    # Any overlap at all?
    if old_base_plasmid_calls & new_base_plasmid_calls:
        return "partial_overlap", old_call.strip()

    # Same replicon family but different variant (e.g. cp32-12 vs cp32-1)
    # These are tie-break differences from identical scores
    old_fam = get_replicon_family(old_call)
    new_fam = get_replicon_family(new_call)
    if old_fam and new_fam and old_fam == new_fam:
        return "same_family_tiebreak", old_call.strip()

    # Genuinely different
    return "different", old_call.strip()


def read_tsv(path: Path, contig_col: str, call_col: str, assembly_col: str = None,
             len_col: str = None):
    """Read a TSV and return dict of (assembly_id, contig_id) -> call.
    If len_col is provided, returns a second dict of (assembly_id, contig_id) -> length."""
    delimiter = "\t" if path.suffix in (".tsv", ".txt") else ","
    records = {}
    lengths = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=delimiter)

        # Validate columns exist
        if contig_col not in reader.fieldnames:
            sys.exit(f"ERROR: Column '{contig_col}' not found in {path.name}. "
                     f"Available: {reader.fieldnames}")
        if call_col not in reader.fieldnames:
            sys.exit(f"ERROR: Column '{call_col}' not found in {path.name}. "
                     f"Available: {reader.fieldnames}")

        has_assembly = assembly_col and assembly_col in reader.fieldnames
        has_len = len_col and len_col in reader.fieldnames

        for row in reader:
            contig = row[contig_col].strip()
            call = row[call_col].strip() if row[call_col] else ""
            assembly = row[assembly_col].strip() if has_assembly else ""
            key = (assembly, contig)
            records[key] = call
            if has_len:
                try:
                    lengths[key] = int(float(row[len_col]))
                except (ValueError, TypeError):
                    lengths[key] = 0

    if len_col:
        return records, lengths
    return records


def _load_overrides(path: Path) -> dict[tuple[str, str], str]:
    """Load manual override file.
    TSV with columns: assembly_id, contig_id, resolved_call"""
    overrides = {}
    delimiter = "\t" if path.suffix in (".tsv", ".txt") else ","
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            assembly = row.get("assembly_id", "").strip()
            contig = row.get("contig_id", "").strip()
            call = row.get("resolved_call", "").strip()
            if assembly and contig:
                overrides[(assembly, contig)] = call
    return overrides


def main():
    ap = argparse.ArgumentParser(
        description="Compare old and new classifier calls and produce homogenized output."
    )
    ap.add_argument("--old", type=Path, required=True,
                    help="Old (published) calls file")
    ap.add_argument("--new", type=Path, required=True,
                    help="New calls file")
    ap.add_argument("--old-contig-col", default="contig_id",
                    help="Contig ID column in old file (default: contig_id)")
    ap.add_argument("--new-contig-col", default="contig_id",
                    help="Contig ID column in new file (default: contig_id)")
    ap.add_argument("--old-call-col", default="final_call",
                    help="Call column in old file (default: final_call)")
    ap.add_argument("--new-call-col", default="final_call",
                    help="Call column in new file (default: final_call)")
    ap.add_argument("--old-assembly-col", default="assembly_id",
                    help="Assembly ID column in old file (default: assembly_id)")
    ap.add_argument("--new-assembly-col", default="assembly_id",
                    help="Assembly ID column in new file (default: assembly_id)")
    ap.add_argument("--output", type=Path, default=Path("comparison.tsv"),
                    help="Output comparison table (default: comparison.tsv)")
    ap.add_argument("--resolved", type=Path, default=None,
                    help="Output resolved/homogenized calls only (optional)")
    ap.add_argument("--review", type=Path, default=None,
                    help="Output detailed review file for 'different' contigs "
                         "with all hits from each database (optional)")
    ap.add_argument("--all-hits-dir", type=Path, default=None,
                    help="Path to classifier output dir containing {db}/tables/*_all.tsv "
                         "(required for --review)")
    ap.add_argument("--min-review-bp", type=int, default=0,
                    help="Only include contigs >= this length in --review output "
                         "(default: 0, no filter)")
    ap.add_argument("--overrides", type=Path, default=None,
                    help="Manual override file (TSV: assembly_id, contig_id, resolved_call). "
                         "These take priority over all other resolution logic.")
    ap.add_argument("--auto-chromosome-bp", type=int, default=0,
                    help="Auto-resolve to new call when new_call is 'chromosome' and "
                         "contig_len >= this threshold (e.g. 100000). 0 = disabled.")
    args = ap.parse_args()

    if args.review and not args.all_hits_dir:
        ap.error("--review requires --all-hits-dir")

    old_calls = read_tsv(args.old, args.old_contig_col, args.old_call_col,
                         args.old_assembly_col)
    new_calls, contig_lengths = read_tsv(args.new, args.new_contig_col, args.new_call_col,
                                          args.new_assembly_col, len_col="contig_len")

    # Load manual overrides
    overrides = {}
    if args.overrides:
        overrides = _load_overrides(args.overrides)
        print(f"Loaded {len(overrides)} manual overrides from {args.overrides}")

    print(f"Old calls: {len(old_calls)} contigs from {args.old.name}")
    print(f"New calls: {len(new_calls)} contigs from {args.new.name}")

    all_keys = sorted(set(old_calls.keys()) | set(new_calls.keys()))

    rows = []
    categories = Counter()
    auto_chrom_count = 0
    override_count = 0

    for key in all_keys:
        assembly, contig = key
        old_call = old_calls.get(key, "")
        new_call = new_calls.get(key, "")
        clen = contig_lengths.get(key, 0)

        if key not in new_calls:
            cat, resolved = "old_only", old_call
        elif key not in old_calls:
            cat, resolved = "new_only", new_call
        else:
            cat, resolved = categorize(old_call, new_call)

        # Auto-chromosome heuristic: trust new call when it says chromosome
        # and contig is large enough ( >100kb)
        # Since chromosomal pf32 paralogs confound calls.
        if (args.auto_chromosome_bp > 0
                and cat == "different"
                and normalize_call(new_call) == "chromosome"
                and clen >= args.auto_chromosome_bp):
            resolved = new_call.strip()
            cat = "auto_chromosome"
            auto_chrom_count += 1

        # Also handle reverse: old says chromosome but new has a specific plasmid call
        # for a small contig — trust the new call
        if (args.auto_chromosome_bp > 0
                and cat == "different"
                and normalize_call(old_call) == "chromosome"
                and not is_empty(new_call)
                and normalize_call(new_call) != "chromosome"
                and clen < args.auto_chromosome_bp):
            resolved = new_call.strip()
            cat = "auto_chromosome"
            auto_chrom_count += 1

        # Manual overrides take final priority
        if key in overrides:
            resolved = overrides[key]
            if cat not in ("exact_match",):
                cat = "manual_override"
                override_count += 1

        categories[cat] += 1
        rows.append({
            "assembly_id": assembly,
            "contig_id": contig,
            "contig_len": clen,
            "old_call": old_call,
            "new_call": new_call,
            "category": cat,
            "resolved_call": resolved,
        })

    # Write comparison table
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t",
                                fieldnames=["assembly_id", "contig_id",
                                            "contig_len", "old_call", "new_call",
                                            "category", "resolved_call"])
        writer.writeheader()
        writer.writerows(rows)

    # Write resolved calls if requested
    if args.resolved:
        with open(args.resolved, "w", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t",
                                    fieldnames=["assembly_id", "contig_id",
                                                "resolved_call"])
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    "assembly_id": row["assembly_id"],
                    "contig_id": row["contig_id"],
                    "resolved_call": row["resolved_call"],
                })

    # Summary
    print(f"\n{'='*50}")
    print(f"{'Category':<25} {'Count':>6}")
    print(f"{'-'*25} {'-'*6}")
    for cat in ["exact_match", "annotation_suffix", "extra_annotation",
                "base_match", "same_family_tiebreak", "auto_chromosome",
                "manual_override", "partial_overlap",
                "different", "new_unclassified", "old_unclassified",
                "old_only", "new_only"]:
        if categories[cat]:
            print(f"  {cat:<23} {categories[cat]:>6}")
    print(f"{'-'*25} {'-'*6}")
    print(f"  {'TOTAL':<23} {sum(categories.values()):>6}")

    if auto_chrom_count:
        print(f"\n  ℹ  {auto_chrom_count} contigs auto-resolved via chromosome heuristic "
              f"(>={args.auto_chromosome_bp}bp)")
    if override_count:
        print(f"  ℹ  {override_count} contigs resolved via manual overrides")

    # Flag items needing manual review
    needs_review = [r for r in rows if r["category"] in
                    ("different", "partial_overlap", "new_unclassified")]
    if needs_review:
        print(f"\n⚠  {len(needs_review)} contigs need manual review:")
        for r in needs_review[:20]:
            print(f"  {r['assembly_id']}/{r['contig_id']} ({r['contig_len']}bp): "
                  f"'{r['old_call']}' -> '{r['new_call']}' [{r['category']}]")
        if len(needs_review) > 20:
            print(f"  ... and {len(needs_review) - 20} more (see {args.output.name})")

    # Detailed review file with all hits
    if args.review and needs_review:
        if args.min_review_bp > 0:
            review_contigs = [r for r in needs_review
                              if r["contig_len"] >= args.min_review_bp]
            skipped = len(needs_review) - len(review_contigs)
            print(f"  Review filter: {skipped} contigs below {args.min_review_bp}bp skipped")
        else:
            review_contigs = needs_review
        if review_contigs:
            _write_review(review_contigs, args.all_hits_dir, args.review)
            print(f"Wrote detailed review ({len(review_contigs)} contigs) -> {args.review}")
        else:
            print(f"  No contigs passed the --min-review-bp filter")

    print(f"\nWrote comparison -> {args.output}")
    if args.resolved:
        print(f"Wrote resolved calls -> {args.resolved}")


def _discover_db_dirs(all_hits_dir: Path) -> dict[str, Path]:
    """Find database subdirectories containing tables/*_all.tsv."""
    dbs = {}
    for child in sorted(all_hits_dir.iterdir()):
        tables_dir = child / "tables"
        if child.is_dir() and tables_dir.is_dir():
            if list(tables_dir.glob("*_all.tsv")):
                dbs[child.name] = tables_dir
    return dbs


def _load_all_hits_for_assembly(tables_dir: Path, assembly_id: str) -> list[str]:
    """Load the all-hits TSV for a given assembly. Returns list of lines."""
    # Try exact match first
    tsv = tables_dir / f"{assembly_id}_all.tsv"
    if tsv.exists():
        with open(tsv) as f:
            return f.readlines()
    # Fallback: glob for partial match
    candidates = list(tables_dir.glob(f"{assembly_id}*_all.tsv"))
    if candidates:
        with open(candidates[0]) as f:
            return f.readlines()
    return []


def _grep_contig(lines: list[str], contig_id: str) -> list[str]:
    """Fixed-string search for contig_id in lines.
    Filters out placeholder rows where hit fields are empty."""
    matches = []
    for l in lines:
        if contig_id not in l:
            continue
        # Check that this is a real hit, not a no-hit placeholder.
        # Placeholder rows have contig_id and contig_len but empty hit fields.
        # Split on tab; a real hit will have plasmid_id populated (column index 3).
        fields = l.rstrip().split("\t")
        if len(fields) > 3 and fields[3].strip():
            matches.append(l)
    return matches


def _write_review(needs_review: list[dict], all_hits_dir: Path, review_path: Path):
    """Write a detailed review file with all hits for each flagged contig."""
    db_dirs = _discover_db_dirs(all_hits_dir)

    if not db_dirs:
        print(f"  WARNING: No database tables found in {all_hits_dir}")
        return

    # Cache loaded files to avoid re-reading
    _cache: dict[tuple[str, str], list[str]] = {}

    with open(review_path, "w") as f:
        f.write(f"# Detailed review of {len(needs_review)} contigs needing manual review\n")
        f.write(f"# Databases found: {', '.join(db_dirs.keys())}\n")
        f.write(f"# Format: all hits per database for each flagged contig\n")
        f.write(f"#\n\n")

        for r in needs_review:
            assembly = r["assembly_id"]
            contig = r["contig_id"]
            old_call = r["old_call"]
            new_call = r["new_call"]
            category = r["category"]

            f.write(f"{'='*100}\n")
            f.write(f"### {assembly} / {contig}\n")
            f.write(f"### old={old_call}  new={new_call}  category={category}\n")
            f.write(f"{'='*100}\n\n")

            for db_name, tables_dir in db_dirs.items():
                cache_key = (db_name, assembly)
                if cache_key not in _cache:
                    _cache[cache_key] = _load_all_hits_for_assembly(tables_dir, assembly)

                lines = _cache[cache_key]
                hits = _grep_contig(lines, contig)

                f.write(f"--- {db_name} ({len(hits)} hits) ---\n")
                if hits:
                    # Write header from first line of file if available
                    if lines and lines[0].startswith("assembly_id"):
                        f.write(f"  {lines[0].rstrip()}\n")
                    for hit in hits:
                        f.write(f"  {hit.rstrip()}\n")
                else:
                    f.write(f"  (no hits)\n")
                f.write(f"\n")

            f.write(f"\n")


if __name__ == "__main__":
    main()
