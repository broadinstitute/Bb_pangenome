#!/usr/bin/env python3
"""
Validate metadb query by comparing Roary output with existing results.

Runs optimized query on Roary groups_to_IDs and compares with
existing clustered_proteins_db_results_pg_v8_sp__v1.2.pkl
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
    # Paths
    database_path = 'metadb/Bbss_db_v3.1.db'
    roary_groups_pkl = 'ref/asm_db/groups_to_IDs_roary_v8_filtered_sp.pkl'
    existing_output = 'metadb/roary_output/clustered_proteins_db_results_pg_v8_sp__v1.2.pkl'
    new_output = 'metadb/roary_output/clustered_proteins_db_results_roary_v3.1_validation.pkl'

    # Load existing output for comparison
    print("Loading existing Roary output for comparison...")
    with open(existing_output, 'rb') as f:
        existing_dict = pickle.load(f)

    print(f"Existing output: {len(existing_dict)} groups")
    existing_non_empty = sum(1 for df in existing_dict.values() if not df.empty)
    existing_total_rows = sum(len(df) for df in existing_dict.values())
    print(f"  Non-empty groups: {existing_non_empty}")
    print(f"  Total rows: {existing_total_rows:,}")

    # Load Roary groups
    print(f"\nLoading Roary groups from: {roary_groups_pkl}")
    with open(roary_groups_pkl, 'rb') as f:
        groups_dict = pickle.load(f)

    total_genes = sum(len(ids) for ids in groups_dict.values())
    print(f"Groups: {len(groups_dict)}, Total genes: {total_genes:,}")

    # Run query with v3.1 database
    print(f"\nQuerying with database: {database_path}")
    new_dict = process_groups(database_path, groups_dict)

    # Save new output
    with open(new_output, 'wb') as f:
        pickle.dump(new_dict, f)
    print(f"\nSaved: {new_output}")

    # Compare results
    new_non_empty = sum(1 for df in new_dict.values() if not df.empty)
    new_total_rows = sum(len(df) for df in new_dict.values())

    print("\n=== COMPARISON ===")
    print(f"{'Metric':<25} {'Existing (v1.2)':<20} {'New (v3.1)':<20}")
    print("-" * 65)
    print(f"{'Total groups':<25} {len(existing_dict):<20} {len(new_dict):<20}")
    print(f"{'Non-empty groups':<25} {existing_non_empty:<20} {new_non_empty:<20}")
    print(f"{'Total annotation rows':<25} {existing_total_rows:<20,} {new_total_rows:<20,}")

    # Check a few specific groups
    print("\n=== SAMPLE GROUP COMPARISON ===")
    sample_groups = list(set(existing_dict.keys()) & set(new_dict.keys()))[:5]
    for group in sample_groups:
        old_rows = len(existing_dict[group])
        new_rows = len(new_dict[group])
        match = "✓" if old_rows == new_rows else "✗"
        print(f"  {group[:30]:<30} old={old_rows:<5} new={new_rows:<5} {match}")
