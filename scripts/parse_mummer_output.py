#!/usr/bin/env python3
import sys
from collections import defaultdict
import pathlib
import pandas as pd
from intervaltree import IntervalTree
from typing import List, Tuple

def simplify_intervals(intervals: List[Tuple[int, int]], sequence_length: int = None) -> int:
    """
    Convert a list of intervals to an interval tree, merge overlaps,
    and return total non-overlapping length.
    
    Args:
        intervals: List of tuples containing (start, end) coordinates
                  Note: Assumes 1-based inclusive coordinates as per MUMmer format
    
    Returns:
        int: Total non-overlapping length of all intervals
    """
    if not intervals:
        return 0
        
    # Convert to interval tree
    # Note: IntervalTree uses half-open intervals [start, end)
    # but MUMmer uses inclusive [start, end], so we add 1 to end
    tree = IntervalTree.from_tuples((start, end + 1, None) for start, end in intervals)
    
    # Merge overlapping intervals
    tree.merge_overlaps()
    
    # Calculate total non-overlapping length
    # Subtract 1 from each interval to convert back to inclusive coordinates
    total_length = sum((interval.end - interval.begin) for interval in sorted(tree))
    
    # If sequence_length is provided, cap the total length
    if sequence_length is not None:
        total_length = min(total_length, sequence_length)
    
    return total_length

def parse_mummer_alignments(coords_file, min_length=100):
    """Parse MUMmer alignments into a network with one edge per sequence pair."""
    # Store alignments with consistent sequence pair ordering
    pair_alignments = defaultdict(lambda: {
        'intervals': [],
        'identity': 0.0,
        'count': 0,
        'seq1_length': 0,
        'seq2_length': 0
    })
    
    with open(coords_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            
            # Skip header or malformed lines
            if len(fields) != 21:
                continue
            try:
                query_id = fields[0]
                ref_id = fields[5]
                
                # Always store with consistent ordering
                seq1, seq2 = sorted([query_id, ref_id])
                
                if seq1 == seq2:  # Skip self alignments
                    continue
                
                # Get coordinates and ensure proper ordering
                qstart, qend = sorted([int(fields[6]), int(fields[7])])
                alignment_length = int(fields[12])
                identity = float(fields[10])
                query_length = int(fields[2])
                ref_length = int(fields[18])
                
                if alignment_length < min_length:
                    continue
                
                # Store data with consistent sequence ordering
                pair_key = (seq1, seq2)
                pair_data = pair_alignments[pair_key]
                
                # Store interval based on which sequence is query
                if query_id == seq1:
                    interval = (qstart, qend)
                    pair_data['seq1_length'] = query_length
                    pair_data['seq2_length'] = ref_length
                else:
                    interval = (qstart, qend)
                    pair_data['seq1_length'] = ref_length
                    pair_data['seq2_length'] = query_length
                
                pair_data['intervals'].append(interval)
                pair_data['identity'] = max(pair_data['identity'], identity)
                pair_data['count'] += 1
                
            except (IndexError, ValueError) as e:
                print(f"Error parsing line: {e}", file=sys.stderr)
                continue
    
    # Create edges from the processed alignments
    edges = []
    for (seq1, seq2), data in pair_alignments.items():
        # Calculate total non-overlapping alignment length
        total_length = simplify_intervals(data['intervals'])
        
        # Ensure we don't exceed sequence lengths
        min_seq_length = min(data['seq1_length'], data['seq2_length'])
        total_length = min(total_length, min_seq_length)
        
        if total_length < min_length:
            continue
        
        edges.append({
            'source': seq1,
            'target': seq2,
            'length': total_length,
            'identity': data['identity'],
            'alignment_count': data['count']
        })
    
    return pd.DataFrame(edges)

def filter_edges(df, min_length=100, min_identity=85):
    """Filter edges based on minimum length and identity."""
    return df[
        (df['length'] >= min_length) & 
        (df['identity'] >= min_identity)
    ]

def print_stats(df):
    """Print network statistics."""
    print("\nNetwork Statistics:")
    nodes = set(df['source']).union(set(df['target']))
    
    print(f"Total edges: {len(df)}")
    print(f"Unique nodes: {len(nodes)}")
    print(f"Average alignments per edge: {df['alignment_count'].mean():.2f}")
    
    print(f"\nAlignment Length Statistics:")
    print(f"Total alignment length: {df['length'].sum():,}")
    print(f"Mean alignment length: {df['length'].mean():.2f}")
    print(f"Median alignment length: {df['length'].median():.2f}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py coords_file output_dir")
        sys.exit(1)
    
    coords_file = pathlib.Path(sys.argv[1])
    output_dir = pathlib.Path(sys.argv[2])
    prefix = coords_file.stem.split('.')[0]
    # Parse and filter alignments
    edges_df = parse_mummer_alignments(coords_file)
    filtered_df = filter_edges(edges_df)
    
    # Write output
    filtered_df.to_csv(output_dir.joinpath(f'{prefix}_homology_network.tsv'), sep='\t', index=False)
    print_stats(filtered_df)

if __name__ == '__main__':
    main()