#!/usr/bin/env python3
"""
Generate B31 gene mappings for Roary v8 output (for comparison with Panaroo).

Roary FASTA headers: >locus_tag gene_name (e.g., >CFHKPG_00001 Uncharacterized protein BB_0001)

Creates:
  - roary_group_to_ID_mapping.txt
  - roary_groupID_to_BBgene.tsv
  - roary_missing_groups.tsv
"""

import os
import pandas as pd

# Paths (relative to Bb_pangenome root)
pg_ref = 'output/results/v8/v8_filtered_sp/roary_v8_filtered_sp/pan_genome_reference.fa'
blast_results = 'output/results/v8/v8_filtered_sp/alignments/pangenome_vs_B31.topHit.tsv'
rs2bb = 'group2BB/B31_refseq_id_to_old_name.tsv'

# Output paths (write to panaroo dir for comparison)
g2idmap = 'panaroo/output/panaroo_out/alignments/roary_group_to_ID_mapping.txt'
g2BBmap = 'panaroo/output/panaroo_out/alignments/roary_groupID_to_BBgene.tsv'
mg_out = 'panaroo/output/panaroo_out/alignments/roary_missing_groups.tsv'

# Step 3: Create group_to_ID_mapping from pan_genome_reference.fa
# For Roary: header is "locus_tag gene_name"
print("Step 3: Creating group_to_ID_mapping (Roary format)...")
with open(pg_ref) as f, open(g2idmap, 'w') as out:
    for line in f:
        if line.startswith('>'):
            line = line.replace('>', '').split(' ')
            seqID = line[0]
            groupID = ' '.join(line[1:]).strip()
            if 'BB' in groupID:
                groupID = 'BB'.join(groupID.split('BB'))
            out.write(f'{seqID}\t{groupID}\n')

with open(g2idmap) as f:
    id_count = sum(1 for _ in f)
print(f"Created: {g2idmap}")
print(f"Gene entries: {id_count}")

# Step 4: Create groupID_to_BBgene mapping
print("\nStep 4: Creating groupID_to_BBgene mapping...")

# Load refseq to BB name mapping
rs2bb_df = pd.read_csv(rs2bb, sep='\t')
rs2bb_df = rs2bb_df.rename(columns={"refseq_id": "gene", "old_name": "old_locus_tag"})
rs2bb_df['old_locus_tag'] = rs2bb_df['old_locus_tag'].apply(lambda x: x.split('%2C')[0])

# Load BLAST results
results = pd.read_csv(blast_results, sep='\t', header=None)
print(f"BLAST results: {len(results)} hits")

# Load ID-to-group mapping
id2group = pd.read_csv(g2idmap, sep='\t', header=None, names=["ID", "group"])
print(f"ID-to-group mappings: {len(id2group)}")

# Create hits DataFrame
hits = pd.DataFrame({
    "ID": results[0],
    "gene": results[1],
    "percent_ident": results[2],
    "alignment_length": results[3],
    "E-score": results[10]
})
hits['gene'] = hits['gene'].apply(lambda x: x.strip('gene-'))

# Merge with group mapping
hits = hits.merge(id2group, on="ID", how="left")
hits = hits.drop_duplicates()
hits = hits.merge(rs2bb_df, on="gene", how="left")

# Get top hits per ID
top_hits = hits.sort_values(
    by=["ID", "group", "gene", "old_locus_tag", "E-score", "percent_ident", "alignment_length"],
    ascending=[True, True, True, False, False, False, False]
).groupby("ID", as_index=False).first()

# Handle duplicate columns
if "group_y" in top_hits.columns:
    top_hits = top_hits.drop(columns=["group_y"]).rename(columns={"group_x": "group"})

# Reorder columns
translatedhits = top_hits[["ID", "group", "gene", "old_locus_tag", "percent_ident", "alignment_length", "E-score"]]

# Find missing groups
missing_groups = id2group[~id2group["group"].isin(translatedhits["group"])]

# Save outputs
translatedhits.to_csv(g2BBmap, sep='\t', header=True, index=False)
missing_groups.to_csv(mg_out, sep='\t', header=True, index=False)

print(f"Created: {g2BBmap}")
print(f"Created: {mg_out}")
print(f"\nSummary:")
print(f"  Total groups: {len(id2group)}")
print(f"  Groups with B31 hits: {len(translatedhits)}")
print(f"  Missing groups: {len(missing_groups)}")
