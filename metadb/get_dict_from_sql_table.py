import sqlite3
import pickle

def read_table_to_dict(db_name, table_name):
    """
    Reads an SQLite table and returns a list of dictionaries containing the data.
    Each dictionary corresponds to one row in the table.
    
    Args:
        db_name (str): Name of the SQLite database file
        table_name (str): Name of the table to read from
        
    Returns:
        list: List of dictionaries containing the gene data
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    # Get all rows from the table
    cursor.execute(f'SELECT * FROM {table_name}')
    
    # Get column names from cursor description
    columns = [description[0] for description in cursor.description]
    
    # Fetch all rows
    rows = cursor.fetchall()
    
    # Convert rows to list of dictionaries
    genes = []
    for row in rows:
        gene_dict = {}
        for column, value in zip(columns, row):
            gene_dict[column] = value
        genes.append(gene_dict)
    
    conn.close()
    return genes

def verify_data(original_dict, retrieved_dict):
    """
    Helper function to verify that all keys from the original dictionary structure
    are present in the retrieved dictionary.
    
    Args:
        original_dict (dict): Original dictionary structure
        retrieved_dict (dict): Retrieved dictionary from database
        
    Returns:
        bool: True if all keys are present, False otherwise
    """
    expected_keys = [
        "feature_type", "gene", "locus_tag", "note", "protein_id",
        "product", "sequence", "replicon", "replicon_name", "start",
        "end", "strand", "db_xrefs", "assembly", "translation",
        "inference", "transl_table"
    ]
    
    return all(key in retrieved_dict for key in expected_keys)
    
def read_and_pickle_genes(db_name: str, table_name: str, pickle_file: str) -> None:
    """
    Reads genes from SQLite database and saves them to a pickle file
    
    Args:
        db_name (str): Name of the SQLite database file
        table_name (str): Name of the table to read from
        pickle_file (str): Name of the pickle file to save to
    """
    # Read from database using our previous function
    genes = read_table_to_dict(db_name, table_name)
    
    # Save to pickle file
    with open(pickle_file, 'wb') as f:
        pickle.dump(genes, f)
    
    return genes

# Example usage:
if __name__ == "__main__":
    db_name = "Bbss_db_v2.db"
    table_name = "ncbi_gb"
    pickle_file = "bakta_genbank_dict_v1.pkl"
    
    # Read from database and save to pickle
    genes = read_and_pickle_genes(db_name, table_name, pickle_file)
    
    # Print some stats to verify
    print(f'Total number of genes retrieved: {len(genes)}')
    
    # Verify first entry has all expected fields
    if genes:
        print("\nFirst gene entry fields:")
        for key in genes[0].keys():
            print(f"- {key}")