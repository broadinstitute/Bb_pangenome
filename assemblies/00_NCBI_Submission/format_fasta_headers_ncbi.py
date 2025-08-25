#!/usr/bin/env python3
"""
Rewrite FASTA headers to NCBI batch-submission format using an identities table.

- Input: directory of FASTA files, plus a single identities.csv describing all assemblies.
- Output: rewritten FASTA files in the output directory.
- Also filters contigs shorter than --min-len (default 200 bp) and logs them.

Usage:
    python rewrite_fasta_ncbi.py input_dir identities.csv output_dir [--min-len 200]

Notes:
- Strain normalization: strip trailing H or P only (e.g., XYZ459H -> XYZ459).
- Chromosome sequences: include [location=chromosome].
- Plasmids: include [plasmid-name=...]; use unique unnamedN when unknown.
- Topology:
    * Use explicit [topology=...] in contig_id if present.
    * Else infer: cp* -> circular; lp*/chromosome -> linear.
- Completeness:
    * If ref_length >= 5000, complete if contig_len within ±7%.
    * Else use name-based ranges (fusion-aware for names like cp32-1+5).
    * Emit [completeness=complete] only when circular and complete.
- Filtering:
    * Drop sequences with length < --min-len and log them to filtered_out.tsv.
"""

from pathlib import Path
import sys
import re
import argparse
import pandas as pd
from Bio import SeqIO

# ---- configuration ----
ORG = "Borreliella burgdorferi"
REF_TOL = 0.07          # ±7% window when a substantial ref_length is available
FULL_REF_MIN = 5000     # consider ref_length "full" if >= this many bp
FUSION_PAD_FRAC = 0.05  # +5% pad on upper bound for fusions
DEFAULT_MIN_LEN = 200   # common NCBI minimum for nucleotide sequences

HDR_RE = re.compile(r'^([A-Za-z0-9]+)_(.+?)_contig_([0-9]+)$')
FUSION_SPLIT = re.compile(r'\s*\+\s*')

def normalize_strain(s: str) -> str:
    """Drop a trailing H or P only (e.g., XYZ459H -> XYZ459; B500P -> B500)."""
    return re.sub(r'(H|P)$', '', s)

def parse_header_token(tok: str):
    """Return (strain_raw, replicon_name, contig_num) from a FASTA id token or None."""
    m = HDR_RE.match(tok)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)

def contig_number_from_id(contig_id_field: str) -> str | None:
    m = re.search(r'contig_(\d+)', contig_id_field)
    return m.group(1) if m else None

def parse_topology_from_contig_id(contig_id_field: str) -> str | None:
    m = re.search(r'\[topology=(\w+)\]', contig_id_field)
    return m.group(1) if m else None

def topology_from_name(name: str) -> str | None:
    if not isinstance(name, str):
        return None
    n = name.lower()
    if n.startswith("cp"):
        return "circular"
    if n.startswith("lp") or n == "chromosome":
        return "linear"
    return None

def base_range(name_lower: str):
    n = name_lower
    if n == "chromosome":    return None
    if n == "cp9":           return (8_000, 9_200)
    if n == "cp26":          return (25_000, 27_500)
    if n.startswith("cp32"): return (28_000, 33_000)
    if n == "lp17":          return (14_000, 21_000)
    if n == "lp25":          return (22_000, 27_000)
    if n.startswith("lp28"): return (27_000, 32_000)
    if n == "lp36":          return (30_000, 39_000)
    if n == "lp38":          return (35_000, 41_000)
    if n == "lp54":          return (50_000, 58_000)
    return None

def name_based_range(plasmid_name: str):
    """Return (lo, hi) bp for a name. Fusion names (e.g., cp32-1+5) sum component ranges and add a small pad."""
    if not isinstance(plasmid_name, str) or not plasmid_name:
        return None
    n = plasmid_name.lower().strip()
    if '+' in n:
        parts = FUSION_SPLIT.split(n)
        lo_sum, hi_sum = 0, 0
        for p in parts:
            r = base_range(p)
            if r is None:
                return None
            lo_sum += r[0]
            hi_sum += r[1]
        hi_sum = int(hi_sum * (1 + FUSION_PAD_FRAC))
        return (lo_sum, hi_sum)
    return base_range(n)

def compute_completeness(contig_len, ref_len, plasmid_name):
    """Return True/False/pd.NA completeness."""
    if pd.isna(contig_len):
        return pd.NA
    # Rule 1: direct ref-length check
    if pd.notna(ref_len) and ref_len >= FULL_REF_MIN:
        lo = ref_len * (1 - REF_TOL)
        hi = ref_len * (1 + REF_TOL)
        if lo <= contig_len <= hi:
            return True
    # Rule 2: fusion-aware name-based range
    rng = name_based_range(plasmid_name)
    if rng:
        lo, hi = rng
        return bool(lo <= contig_len <= hi)
    return pd.NA

def build_index_from_table(tcsv: Path):
    """Build lookup keyed by (assembly_id, contig_num)."""
    df = pd.read_csv(tcsv)
    df["contig_num"] = df["contig_id"].map(contig_number_from_id)
    df["topology_explicit"] = df["contig_id"].map(parse_topology_from_contig_id)
    df["topology_inferred"] = df["plasmid_name"].map(topology_from_name)
    df["complete_flag"] = df.apply(
        lambda r: compute_completeness(r["contig_len"], r["ref_length"], r["plasmid_name"]), axis=1
    )

    idx = {}
    for r in df.itertuples(index=False):
        key = (r.assembly_id, str(r.contig_num))
        idx[key] = dict(
            plasmid_name=r.plasmid_name,
            topology_explicit=r.topology_explicit,
            topology_inferred=r.topology_inferred,
            complete=r.complete_flag,
        )
    return idx

def rewrite_and_filter_fasta(fasta_in: Path, fasta_out: Path, table_idx: dict,
                             min_len: int, filtered_log_rows: list[dict]):
    """Rewrite headers for one FASTA file; filter short contigs and log them."""
    unnamed_counter: dict[str, int] = {}

    def next_unnamed(asm: str) -> str:
        unnamed_counter[asm] = unnamed_counter.get(asm, 0) + 1
        return f"unnamed{unnamed_counter[asm]}"

    written = 0
    with fasta_out.open("w") as out_handle:
        for rec in SeqIO.parse(str(fasta_in), "fasta"):
            # Filter by contig length
            rec_len = len(rec.seq)
            if rec_len < min_len:
                # Record in filtered log
                filtered_log_rows.append({
                    "source_file": fasta_in.name,
                    "original_id": rec.id,
                    "length": rec_len,
                    "reason": f"< min_len ({min_len})"
                })
                continue

            parsed = parse_header_token(rec.id)
            if not parsed:
                # Pass-through if header doesn't match expected pattern
                SeqIO.write(rec, out_handle, "fasta")
                written += 1
                continue

            strain_raw, rep_name, contig_num = parsed
            strain_out = normalize_strain(strain_raw)
            assembly_id = strain_out  # aligns with identities.csv (e.g., B331, B500)

            meta = table_idx.get((assembly_id, contig_num), {})
            plasmid_name = meta.get("plasmid_name", rep_name)

            # Unknown plasmids -> unnamedN per assembly
            p_lower = str(plasmid_name).lower()
            if p_lower in ("unknownreplicon", "unknown", "nan", ""):
                plasmid_name = next_unnamed(assembly_id)

            topo = meta.get("topology_explicit") or meta.get("topology_inferred") or topology_from_name(plasmid_name)
            complete = meta.get("complete")

            # Build NCBI-compliant header
            new_id = f"contig_{contig_num}"
            desc_bits = [f"[organism={ORG}]", f"[isolate={strain_out}]"]

            if str(plasmid_name).lower() == "chromosome":
                desc_bits.append("[location=chromosome]")
                if topo in ("linear", "circular"):
                    desc_bits.append(f"[topology={topo}]")
                if complete is True: 
                    desc_bits.append("[completeness=complete]")
            else:
                desc_bits.append(f"[plasmid-name={plasmid_name}]")


            rec.id = new_id
            rec.name = new_id
            rec.description = " ".join(desc_bits)
            SeqIO.write(rec, out_handle, "fasta")
            written += 1
    return written

def find_fastas(input_dir: Path):
    """Yield FASTA-like files. Adjust patterns as needed."""
    exts = (".fa", ".fna", ".fasta", ".fas")
    for p in sorted(input_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in exts:
            yield p

def main():
    ap = argparse.ArgumentParser(description="Rewrite FASTA headers to NCBI format; filter short contigs.")
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("identities_csv", type=Path)
    ap.add_argument("output_dir", type=Path)
    ap.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN,
                    help=f"Minimum contig length to keep (default: {DEFAULT_MIN_LEN})")
    args = ap.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"Input directory not found: {args.input_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    table_idx = build_index_from_table(args.identities_csv)

    filtered_log_rows: list[dict] = []
    total_in, total_out, total_filtered = 0, 0, 0

    for fasta_in in find_fastas(args.input_dir):
        fasta_out = args.output_dir / fasta_in.name
        written = rewrite_and_filter_fasta(
            fasta_in, fasta_out, table_idx, args.min_len, filtered_log_rows
        )
        # count input records
        in_count = sum(1 for _ in SeqIO.parse(str(fasta_in), "fasta"))
        total_in += in_count
        total_out += written
        total_filtered += (in_count - written)
        print(f"[OK] {fasta_in.name}: kept {written}/{in_count} (filtered {in_count - written}) -> {fasta_out.name}")

    # write filtered-out log
    if filtered_log_rows:
        log_path = args.output_dir / "filtered_out.tsv"
        pd.DataFrame(filtered_log_rows).to_csv(log_path, sep="\t", index=False)
        print(f"[LOG] Filtered contigs written to {log_path}")

    print(f"\nSummary: kept {total_out}/{total_in} total (filtered {total_filtered}, min_len={args.min_len})")

if __name__ == "__main__":
    main()

