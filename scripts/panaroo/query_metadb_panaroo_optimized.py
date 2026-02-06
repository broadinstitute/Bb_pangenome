#!/usr/bin/env python3
"""
Query metadb for Panaroo clustered proteins results (optimized version).

Optimizations over original:
  - Single DB connection (reused)
  - Batch queries using IN clause instead of per-gene queries
  - Chunked processing to avoid SQL query length limits

Creates:
  - clustered_proteins_db_results_panaroo.pkl (Step 6)
"""

import pandas as pd
import pickle
import sqlite3
from tqdm import tqdm
import os

def query_genes_batch(conn, locus_tags, chunk_size=500):
    """Query multiple locus_tags in batches using IN clause."""
    all_results = []

    for i in range(0, len(locus_tags), chunk_size):
        chunk = locus_tags[i:i + chunk_size]
        placeholders = ','.join(['?' for _ in chunk])
        query = f"""
            SELECT * FROM annotations
            WHERE feature_type = 'CDS'
            AND locus_tag IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=chunk)
        all_results.append(df)

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


def process_groups(db_path, groups_dict):
    """Process all gene groups with a single DB connection."""
    conn = sqlite3.connect(db_path)

    result_dict = {}

    for group, gene_ids in tqdm(groups_dict.items(), desc="Processing groups"):
        if not gene_ids:
            result_dict[group] = pd.DataFrame()
            continue

        df = query_genes_batch(conn, gene_ids)
        result_dict[group] = df

    conn.close()
    return result_dict


if __name__ == "__main__":
    # Paths (relative to Bb_pangenome root)
    database_path = 'metadb/Bbss_db_v3.1.db'
    groups_pkl = 'panaroo/output/panaroo_out/groups_to_IDs_panaroo.pkl'
    output_pkl = 'panaroo/output/panaroo_out/clustered_proteins_db_results_panaroo.pkl'

    # Verify database exists
    if not os.path.exists(database_path):
        print(f"ERROR: Database not found at {database_path}")
        exit(1)

    print(f"Database: {database_path}")
    print(f"Loading groups from: {groups_pkl}")

    with open(groups_pkl, 'rb') as jar:
        groups_dict = pickle.load(jar)

    total_genes = sum(len(ids) for ids in groups_dict.values())
    print(f"Processing {len(groups_dict)} gene groups ({total_genes:,} total genes)...")

    # Check for index on locus_tag (informational)
    conn = sqlite3.connect(database_path)
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='index'")
    indexes = [r[0] for r in cur.fetchall()]
    conn.close()

    if not any('locus' in idx.lower() for idx in indexes):
        print("Note: No index found on locus_tag - consider adding one for faster queries")
        print("  CREATE INDEX idx_locus_tag ON annotations(locus_tag);")

    group_dfs_dict = process_groups(database_path, groups_dict)

    with open(output_pkl, 'wb') as jar:
        pickle.dump(group_dfs_dict, jar)

    print(f"\nCreated: {output_pkl}")

    non_empty = sum(1 for df in group_dfs_dict.values() if not df.empty)
    total_rows = sum(len(df) for df in group_dfs_dict.values())
    print(f"Groups with results: {non_empty}/{len(group_dfs_dict)}")
    print(f"Total annotation rows: {total_rows:,}")
