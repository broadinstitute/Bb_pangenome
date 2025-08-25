#!/usr/bin/env python3
import sys
from Bio import SeqIO
import numpy as np
import pandas as pd

def get_contig_lengths(fasta_file):
    """Extract contig lengths from FASTA file."""
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths

def create_alignment_matrix(network_file, fasta_file):
    """Create alignment length matrix from network file and FASTA file."""
    # Read network data
    network = pd.read_csv(network_file, sep='\t')
    
    # Get contig lengths from FASTA
    contig_lengths = get_contig_lengths(fasta_file)
    
    # Get unique contigs from both sources
    all_contigs = sorted(set(network['source']).union(set(network['target'])))
    n_contigs = len(all_contigs)
    
    # Create contig to index mapping
    contig_to_idx = {contig: i for i, contig in enumerate(all_contigs)}
    
    # Initialize matrix with zeros
    matrix = np.zeros((n_contigs, n_contigs))
    
    # Fill diagonal with actual contig lengths from FASTA
    for contig, length in contig_lengths.items():
        if contig in contig_to_idx:
            idx = contig_to_idx[contig]
            matrix[idx, idx] = length
    
    # Fill off-diagonal elements with alignment lengths
    for _, row in network.iterrows():
        i = contig_to_idx[row['source']]
        j = contig_to_idx[row['target']]
        matrix[i, j] = row['length']
        matrix[j, i] = row['length']  # Make symmetric
    
    matrix_df = pd.DataFrame(matrix, index=all_contigs, columns=all_contigs)
    return matrix_df, contig_lengths

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py homology_network.tsv input.fasta")
        sys.exit(1)
    
    network_file = sys.argv[1]
    fasta_file = sys.argv[2]
    
    # Create and save matrix
    matrix_df, contig_lengths = create_alignment_matrix(network_file, fasta_file)
    matrix_df.to_csv('alignment_matrix.tsv', sep='\t')
    
    # Print matrix statistics
    print("\nMatrix Statistics:")
    print(f"Matrix dimensions: {matrix_df.shape[0]} x {matrix_df.shape[1]}")
    diag_sum = np.trace(matrix_df.values)
    print(f"Sum of diagonal elements (total contig lengths): {diag_sum:,.0f}")
    off_diag = matrix_df.values[~np.eye(matrix_df.shape[0], dtype=bool)]
    print(f"Number of non-zero off-diagonal elements: {np.count_nonzero(off_diag)}")
    print(f"Average non-zero alignment length: {np.mean(off_diag[off_diag > 0]):.2f}")
    
    # Print any contigs in network but missing from FASTA
    network_contigs = set(matrix_df.index)
    fasta_contigs = set(contig_lengths.keys())
    missing = network_contigs - fasta_contigs
    if missing:
        print("\nWarning: Found contigs in network but missing from FASTA:")
        for contig in sorted(missing)[:5]:
            print(f"  {contig}")
        if len(missing) > 5:
            print(f"  ...and {len(missing)-5} more")

if __name__ == '__main__':
    main()