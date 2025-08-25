#!/bin/bash

INPUT="nucl_v5"
OUTPUT_FILE="all_alignments.tsv"
PROC_COUNT=0
clear
echo "processing all alignments and concatenating to a single file."
find $INPUT -name "*.tsv" -type f | while read file; do
    cat $file >> $OUTPUT_FILE 
    echo >> $OUTPUT_FILE
    # use '\033[2K' to overwrite the previous output and keep on one line.
    echo -ne "\rCount: $PROC_COUNT    "
    ((PROC_COUNT+=1))
done
echo "Done!"
exit 0
