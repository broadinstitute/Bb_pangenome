#!/usr/bin/env bash

INPUT=$1
OUTPUT=$2
THREADS="$((nproc - 2))"
mkdir -p "$OUTPUT"

CMD="docker run --rm -v $(pwd):/opt -it mjfos2r/panaroo:latest -c "
CMD+="panaroo -t $THREADS -i \"${INPUT}\"/*.gff3 -o $OUTPUT --clean-mode strict"


# Need to add the sketch dbs to the container....
#QC_CMD="panaroo-qc -t $THREADS --graph_type all -i ${INPUT}/*.gff"
# run QC
#echo "Running panaroo-qc, please stand by..."

echo "Running panaroo. please stand by..."
echo "$CMD"
$CMD

