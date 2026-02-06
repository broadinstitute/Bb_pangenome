#!/usr/bin/env python3
"""
Query metadb for Panaroo clustered proteins results.

Adapted from: metadb/query_database_clustered_proteins.py

Creates:
  - clustered_proteins_db_results_panaroo.pkl (Step 6)
"""

import pandas as pd
import pickle
import sqlite3
from tqdm import tqdm
import multiprocessing
import os

def execute_query(db_path, query):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

def process_gene_group(db_path, group, gene_ids, progress_queue):
    group_df = pd.DataFrame()
    for gene in gene_ids:
        query = f"SELECT * FROM annotations WHERE feature_type = 'CDS' AND locus_tag = '{gene}';"
        result_df = execute_query(db_path, query)
        group_df = pd.concat([group_df, result_df], ignore_index=True)
        progress_queue.put(1)
    return group, group_df

def update_progress_bar(total_tasks, progress_queue):
    with tqdm(total=total_tasks, desc="Querying genes") as pbar:
        for _ in range(total_tasks):
            progress_queue.get()
            pbar.update(1)

def parallel_process_groups(db_path, groups_dict):
    manager = multiprocessing.Manager()
    progress_queue = manager.Queue()
    total_tasks = sum(len(ids) for ids in groups_dict.values())

    print(f"Total genes to query: {total_tasks}")

    progress_process = multiprocessing.Process(
        target=update_progress_bar,
        args=(total_tasks, progress_queue)
    )
    progress_process.start()

    num_workers = max(1, multiprocessing.cpu_count() - 2)
    print(f"Using {num_workers} worker processes")

    with multiprocessing.Pool(processes=num_workers) as pool:
        args = [(db_path, gene_group, ids, progress_queue)
                for gene_group, ids in groups_dict.items()]
        results = pool.starmap(process_gene_group, args)

    result_dict = {gene_group: df for gene_group, df in results}
    progress_process.join()
    return result_dict

if __name__ == "__main__":
    # Database path (relative to Bb_pangenome root)
    database_path = 'metadb/Bbss_db_v3.1.db'

    # Input/output paths (relative to Bb_pangenome)
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

    print(f"Processing {len(groups_dict)} gene groups...")

    group_dfs_dict = parallel_process_groups(database_path, groups_dict)

    with open(output_pkl, 'wb') as jar:
        pickle.dump(group_dfs_dict, jar)

    print(f"\nCreated: {output_pkl}")

    non_empty = sum(1 for df in group_dfs_dict.values() if not df.empty)
    print(f"Groups with results: {non_empty}/{len(group_dfs_dict)}")
