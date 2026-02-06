#!/usr/bin/env python3
"""
Generate all intermediate files from Panaroo output.

Consolidates Steps 2-6 into a single script:
  - BLAST pan_genome_reference.fa against B31
  - Create group_to_ID_mapping
  - Create groupID_to_BBgene mapping
  - Create groups_to_IDs pickle
  - Query metadb for clustered_proteins

Usage:
  uv run python3 panaroo/scripts/generate_all_intermediates.py <panaroo_output_dir>

Example:
  uv run python3 panaroo/scripts/generate_all_intermediates.py panaroo/debugging_2/no_merge_paralogs
"""

import os
import sys
import subprocess
import pickle
import sqlite3
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def run_blast(panaroo_dir: Path, output_dir: Path, blast_db: str = "group2BB/blastDB/B31_prot"):
    """Step 2: BLAST pan_genome_reference against B31."""
    print("\n[Step 2] Running BLAST against B31...")

    pan_ref = panaroo_dir / "pan_genome_reference.fa"
    blast_out = output_dir / "panaroo_vs_B31.topHit.tsv"

    if not pan_ref.exists():
        raise FileNotFoundError(f"pan_genome_reference.fa not found in {panaroo_dir}")

    cmd = [
        "blastn",
        "-task", "blastn",
        "-query", str(pan_ref),
        "-db", blast_db,
        "-outfmt", "6 qseqid sseqid pident length evalue",
        "-max_target_seqs", "1",
        "-max_hsps", "1",
        "-evalue", "1e-10",
        "-num_threads", "8",
        "-out", str(blast_out)
    ]

    subprocess.run(cmd, check=True)

    # Count results
    with open(blast_out) as f:
        n_hits = sum(1 for _ in f)
    print(f"  BLAST complete: {n_hits} hits")

    return blast_out


def create_group_mappings(panaroo_dir: Path, output_dir: Path, blast_results: Path):
    """Steps 3 & 4: Create group_to_ID and groupID_to_BBgene mappings."""
    print("\n[Steps 3-4] Creating group mappings...")

    pan_ref = panaroo_dir / "pan_genome_reference.fa"
    rs2bb_file = Path("group2BB/B31_refseq_id_to_old_name.tsv")

    g2id_out = output_dir / "panaroo_group_to_ID_mapping.txt"
    g2bb_out = output_dir / "panaroo_groupID_to_BBgene.tsv"
    missing_out = output_dir / "panaroo_missing_groups.tsv"

    # Step 3: For Panaroo, gene name IS both ID and group
    with open(pan_ref) as f, open(g2id_out, 'w') as out:
        for line in f:
            if line.startswith('>'):
                gene_name = line[1:].strip()
                out.write(f'{gene_name}\t{gene_name}\n')

    # Load RefSeq to BB mapping
    rs2bb = {}
    with open(rs2bb_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                rs2bb[parts[0]] = parts[1]

    # Load group_to_ID mapping
    group2id = {}
    with open(g2id_out) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                group2id[parts[0]] = parts[1]

    # Load BLAST results
    blast_hits = {}
    with open(blast_results) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                qseqid, sseqid, pident, length, evalue = parts[:5]
                blast_hits[qseqid] = (sseqid, pident, length, evalue)

    # Step 4: Create groupID_to_BBgene mapping
    mapped = 0
    missing = 0

    with open(g2bb_out, 'w') as out, open(missing_out, 'w') as miss:
        out.write("ID\tgroup\tgene\told_locus_tag\tpercent_ident\talignment_length\tE-score\n")
        miss.write("group\treason\n")

        for group, gene_id in group2id.items():
            if gene_id in blast_hits:
                sseqid, pident, length, evalue = blast_hits[gene_id]
                # Handle different formats: gene-BB_RS00010, BB_RS00010, or with pipes
                rs_id = sseqid.split('|')[0] if '|' in sseqid else sseqid
                rs_id = rs_id.replace('gene-', '')  # Strip gene- prefix if present

                if rs_id in rs2bb:
                    old_tag = rs2bb[rs_id]
                    out.write(f"{gene_id}\t{group}\t{rs_id}\t{old_tag}\t{pident}\t{length}\t{evalue}\n")
                    mapped += 1
                else:
                    miss.write(f"{group}\tno_BB_mapping_for_{rs_id}\n")
                    missing += 1
            else:
                miss.write(f"{group}\tno_blast_hit\n")
                missing += 1

    print(f"  Mapped to B31: {mapped}")
    print(f"  Missing: {missing}")

    return g2id_out, g2bb_out


def create_groups_pickle(panaroo_dir: Path, output_dir: Path):
    """Step 5: Create groups_to_IDs pickle from gene_presence_absence.csv."""
    print("\n[Step 5] Creating groups_to_IDs pickle...")

    gpa_csv = panaroo_dir / "gene_presence_absence.csv"
    output_pkl = output_dir / "groups_to_IDs.pkl"

    df = pd.read_csv(gpa_csv)
    sample_cols = df.columns[3:]  # Skip Gene, Non-unique Gene name, Annotation

    groups_to_genelists = {}
    for _, row in df.iterrows():
        group_id = row['Gene']
        gene_list = []
        for col in sample_cols:
            val = row[col]
            if pd.notna(val) and val != '':
                genes = str(val).split(';')
                gene_list.extend([g.strip() for g in genes if g.strip()])
        groups_to_genelists[group_id] = gene_list

    with open(output_pkl, 'wb') as jar:
        pickle.dump(groups_to_genelists, jar)

    total_genes = sum(len(ids) for ids in groups_to_genelists.values())
    print(f"  Groups: {len(groups_to_genelists)}")
    print(f"  Total gene IDs: {total_genes:,}")

    return output_pkl


def query_metadb(groups_pkl: Path, output_dir: Path, db_path: str = "metadb/Bbss_db_v3.1.db"):
    """Step 6: Query metadb for clustered_proteins results."""
    print("\n[Step 6] Querying metadb...")

    output_pkl = output_dir / "clustered_proteins_db_results.pkl"

    with open(groups_pkl, 'rb') as f:
        groups_dict = pickle.load(f)

    conn = sqlite3.connect(db_path)
    result_dict = {}

    for group, gene_ids in tqdm(groups_dict.items(), desc="  Processing groups"):
        if not gene_ids:
            result_dict[group] = pd.DataFrame()
            continue

        # Batch query
        all_results = []
        chunk_size = 500
        for i in range(0, len(gene_ids), chunk_size):
            chunk = gene_ids[i:i + chunk_size]
            placeholders = ','.join(['?' for _ in chunk])
            query = f"""
                SELECT * FROM annotations
                WHERE feature_type = 'CDS'
                AND locus_tag IN ({placeholders})
            """
            df = pd.read_sql_query(query, conn, params=chunk)
            all_results.append(df)

        if all_results:
            result_dict[group] = pd.concat(all_results, ignore_index=True)
        else:
            result_dict[group] = pd.DataFrame()

    conn.close()

    with open(output_pkl, 'wb') as f:
        pickle.dump(result_dict, f)

    non_empty = sum(1 for df in result_dict.values() if not df.empty)
    total_rows = sum(len(df) for df in result_dict.values())
    print(f"  Groups with results: {non_empty}/{len(result_dict)}")
    print(f"  Total annotation rows: {total_rows:,}")

    return output_pkl


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    panaroo_dir = Path(sys.argv[1])

    if not panaroo_dir.exists():
        print(f"Error: Directory not found: {panaroo_dir}")
        sys.exit(1)

    # Create output subdirectory for intermediates
    output_dir = panaroo_dir / "intermediates"
    output_dir.mkdir(exist_ok=True)

    print(f"=" * 60)
    print(f"Generating intermediates for: {panaroo_dir}")
    print(f"Output directory: {output_dir}")
    print(f"=" * 60)

    # Run all steps
    blast_results = run_blast(panaroo_dir, output_dir)
    create_group_mappings(panaroo_dir, output_dir, blast_results)
    groups_pkl = create_groups_pickle(panaroo_dir, output_dir)
    query_metadb(groups_pkl, output_dir)

    print(f"\n" + "=" * 60)
    print(f"All intermediates generated in: {output_dir}")
    print(f"=" * 60)

    # List outputs
    print("\nGenerated files:")
    for f in sorted(output_dir.iterdir()):
        size = f.stat().st_size
        if size > 1_000_000:
            size_str = f"{size / 1_000_000:.1f}M"
        elif size > 1000:
            size_str = f"{size / 1000:.1f}K"
        else:
            size_str = f"{size}B"
        print(f"  {f.name}: {size_str}")


if __name__ == "__main__":
    main()
