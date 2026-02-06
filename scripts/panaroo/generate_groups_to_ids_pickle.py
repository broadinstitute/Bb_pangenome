#!/usr/bin/env python3
"""
Parse Panaroo gene_presence_absence.csv to create groups_to_IDs pickle.

Adapted from: notebooks/metadb/generate_groups_to_geneIDs.ipynb

Creates:
  - groups_to_IDs_panaroo.pkl (Step 5)
"""

import pandas as pd
import pickle

# Paths (relative to Bb_pangenome root)
gpa_csv = 'panaroo/output/panaroo_out/gene_presence_absence.csv'
output_pkl = 'panaroo/output/panaroo_out/groups_to_IDs_panaroo.pkl'

print("Step 5: Creating groups_to_IDs pickle...")

# Read Panaroo gene_presence_absence.csv
df = pd.read_csv(gpa_csv)
print(f"Loaded {len(df)} gene groups from CSV")

# First 3 columns are: Gene, Non-unique Gene name, Annotation
# Remaining columns are sample names with locus_tags as values
sample_cols = df.columns[3:]
print(f"Sample columns: {len(sample_cols)}")

groups_to_genelists = {}
for _, row in df.iterrows():
    group_id = row['Gene']
    gene_list = []
    for col in sample_cols:
        val = row[col]
        if pd.notna(val) and val != '':
            # Handle multiple genes per cell (semicolon-separated in Panaroo)
            genes = str(val).split(';')
            gene_list.extend([g.strip() for g in genes if g.strip()])
    groups_to_genelists[group_id] = gene_list

# Save pickle
with open(output_pkl, 'wb') as jar:
    pickle.dump(groups_to_genelists, jar)

print(f"\nCreated: {output_pkl}")
print(f"Total gene groups: {len(groups_to_genelists)}")

# Stats
total_genes = sum(len(v) for v in groups_to_genelists.values())
print(f"Total gene instances: {total_genes}")
print(f"Average genes per group: {total_genes / len(groups_to_genelists):.1f}")
