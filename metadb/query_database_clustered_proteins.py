#!/usr/local/bin/python3.11

import pandas
import pickle
import sqlite3
from tqdm import tqdm
import multiprocessing

def execute_query(db_path, query):
    """Executes a given SQL query using a database connection.

    Parameters:
        db_path: Path to the SQLite database file.
        query: SQL query to execute.
    Returns:
        DataFrame containing the results of the query.
    """

    conn = sqlite3.connect(db_path)
    df = pandas.read_sql_query(query, conn)
    conn.close()
    return df

def process_gene_group(db_path, group, gene_ids, progress_queue):
    """Processes gene_ids for gene groups output by Roary pangenome pipeline. A group is fed into this function and each gene_id within the list of gene_ids is
    queried against the given SQLite3 database. A concatenated dataframe is returned.

    Parameters:
        db_path: Path to the SQLite database file.
        gene_group: The gene_group name.
        gene_ids: List of gene_ids for the given group.
        progress_queue: Queue for global progress tracking.

    Returns:
        Tuple of group name and concatenated DataFrame for the group.
    """
    group_df = pandas.DataFrame() # set up empty df
    for gene in gene_ids: # using the list of gene_ids for this group, iterate through and query the db.
        query = f"SELECT * FROM annotations WHERE feature_type = 'CDS' AND locus_tag = '{gene}';"  # Adjust the query as needed
        result_df = execute_query(db_path, query) # execute the query, save df as result
        group_df = pandas.concat([group_df, result_df], ignore_index=True) # concat results into single group_df
        progress_queue.put(1) # put a value into the progress queue so we can update our pbar!
    return group, group_df # return group_id and group_df

def update_progress_bar(total_tasks, progress_queue):
    """Monitors the progress queue and updates the progress bar.

    Parameters:
        total_tasks: Total number of tasks to complete.
        progress_queue: Queue to receive progress updates.
    """
    with tqdm(total=total_tasks) as pbar:
        for _ in range(total_tasks):
            progress_queue.get()
            pbar.update(1)

def parallel_process_groups(db_path, groups_dict):
    """Processes all groups in parallel by distributing each group to a separate worker.

    Parameters
        db_path: Path to the SQLite database file.
        groups_dict: Dictionary of gene_groups with their associated list of gene_ids.
    Returns:
        Dictionary of groups with their concatenated DataFrames.
    """
    manager = multiprocessing.Manager() # set up SyncManager "To manage shared state"
    progress_queue = manager.Queue() # instantiate process queue w/in SyncManager
    total_tasks = sum(len(ids) for ids in groups_dict.values()) # Count how much we gotta do.

    # Start the progress monitor process
    progress_process = multiprocessing.Process(target=update_progress_bar, args=(total_tasks, progress_queue))
    progress_process.start()

    with multiprocessing.Pool(processes=(multiprocessing.cpu_count()-2)) as pool:
        # Create a list of arguments for each group
        args = [(db_path, gene_group, ids, progress_queue) for gene_group, ids in groups_dict.items()]
        # Execute the processing in parallel
        results = pool.starmap(process_gene_group, args)

    # Convert the list of results to a dictionary
    result_dict = {gene_group: df for gene_group, df in results}

    # Wait for the progress monitor thread to finish
    progress_process.join()

    return result_dict

if __name__ == "__main__":
    # define I/O
    database_path = '../ref/asm_db/Bbss_db_v3.1.db' # /home/mf019/longread_pangenome/expanded_dataset_analysis/
    groups_to_genelists_pickle = '../ref/asm_db/groups_to_IDs_roary_v8_filtered_nsp.pkl' # /home/mf019/longread_pangenome/expanded_dataset_analysis
    clustered_protein_results_db_pkl = 'clustered_proteins_db_results_pg_v8_nsp__v1.pkl' # /home/mf019/longread_pangenome/expanded_dataset_analysis/metadb/
    # load pickle
    with open(groups_to_genelists_pickle, 'rb') as jar:
        groups_dict = pickle.load(jar)

    # Process the groups in parallel
    group_dfs_dict = parallel_process_groups(database_path, groups_dict)
    with open(clustered_protein_results_db_pkl, 'wb') as jar:
        pickle.dump(group_dfs_dict, jar)

    # Print the results
    for group, df in group_dfs_dict.items():
        print(f"num of results for {group}:")
        print(len(df))