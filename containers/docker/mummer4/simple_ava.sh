#!/bin/bash

# Input FASTA file containing all contigs
INPUT_FASTA=$1
OUTPUT_PREFIX=$2
THREADS=$(($(nproc)-2))

# Run nucmer for self-alignment
# removed --nosimplify \ # as it appears to be creating difficulties downstream...
time nucmer \
  --maxmatch \
  -t "$THREADS" \
  -p "${OUTPUT_PREFIX}" \
  "${INPUT_FASTA}" \
  "${INPUT_FASTA}"

# Convert to coords format showing ALL alignments
show-coords -c -r -l -T "${OUTPUT_PREFIX}.delta" > "${OUTPUT_PREFIX}.coords"

# Tab-delimited format for easier parsing
show-coords -r -B "${OUTPUT_PREFIX}.delta" > "${OUTPUT_PREFIX}.coords.tab"
