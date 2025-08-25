#!/usr/bin/env bash
set -euo pipefail

ACC_FILE="accessions.txt"            # one SRR/ERR/DRR per line
GENOME_SIZE_BP=1400000

# collapse into comma-separated list
ACC_CSV=$(paste -sd, "$ACC_FILE")

# fetch run info from ENA
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC_CSV}&result=read_run&fields=run_accession,sample_accession,base_count,read_count,library_layout&format=tsv" \
| awk -v g="$GENOME_SIZE_BP" 'BEGIN{FS=OFS="\t"}
    NR==1 {print $0,"coverage_x"; next}
    {
        b = ($3==""?0:$3)+0
        cov = (g>0 ? b/g : 0)
        print $0, sprintf("%.2f", cov)
        total += b
    }
    END {
        if (NR>1) print "TOTAL","","",total,"",sprintf("%.2f",total/g)
    }'

