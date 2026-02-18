#!/usr/bin/env bash
#
# 04_batch_submit.sh — Submit all assemblies to ENA production.
#
# Usage:
#   ./04_batch_submit.sh <manifest_dir> <webin_user> <webin_password> [webin_cli_jar]
#
# Submits all manifests to ENA production. Logs everything.
# Collects receipt information for accession tracking.
#
# IMPORTANT: This creates REAL accessions. Run validation and test first!

set -euo pipefail

MANIFEST_DIR="${1:?Usage: $0 <manifest_dir> <webin_user> <webin_password> [webin_cli_jar]}"
WEBIN_USER="${2:?Provide Webin username}"
#why is this plaintext???
#WEBIN_PASS="${3:?Provide Webin password}"
stty -echo
printf "Webin password: "
read WEBIN_PASS
stty echo
echo
WEBIN_CLI="${3:-webin-cli.jar}"
LOG_DIR="${MANIFEST_DIR}/logs/submit"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RECEIPT_LOG="${MANIFEST_DIR}/submission_receipts_${TIMESTAMP}.tsv"

mkdir -p "$LOG_DIR"

# Initialize receipt log
echo -e "assembly_name\tstatus\tlog_file\ttimestamp" > "$RECEIPT_LOG"

TOTAL=0
PASS=0
FAIL=0
SKIP=0

# Count total manifests
for manifest in "$MANIFEST_DIR"/manifest_*.txt; do
    [ -f "$manifest" ] && TOTAL=$((TOTAL + 1))
done

echo "============================================================"
echo "PRODUCTION submission: $TOTAL assemblies"
echo "  Started: $(date)"
echo "  Receipt log: $RECEIPT_LOG"
echo "============================================================"
echo ""

# Confirmation prompt
read -p "This will create REAL accessions. Continue? (yes/no): " confirm
if [ "$confirm" != "yes" ]; then
    echo "Aborted."
    exit 0
fi

COUNT=0
for manifest in "$MANIFEST_DIR"/manifest_*.txt; do
    [ -f "$manifest" ] || continue
    COUNT=$((COUNT + 1))

    sample=$(basename "$manifest" .txt | sed 's/^manifest_//')
    logfile="$LOG_DIR/${sample}.log"

    echo -n "  [$COUNT/$TOTAL] Submitting: $sample ... "
    mkdir -p "$LOG_DIR/$sample"
    if java -jar "$WEBIN_CLI" \
        -context genome \
        -userName "$WEBIN_USER" \
        -password "$WEBIN_PASS" \
        -manifest "$manifest" \
        -outputDir "$LOG_DIR/$sample/" \
        -submit \
        > "$logfile" 2>&1; then
        echo "OK"
        echo -e "${sample}\tSUCCESS\t${logfile}\t$(date -Iseconds)" >> "$RECEIPT_LOG"
        PASS=$((PASS + 1))
    else
        echo "FAIL (see $logfile)"
        echo -e "${sample}\tFAILED\t${logfile}\t$(date -Iseconds)" >> "$RECEIPT_LOG"
        FAIL=$((FAIL + 1))
    fi

    # Be respectful to ENA servers
    sleep 5
done

echo ""
echo "============================================================"
echo "Submission complete: $PASS succeeded, $FAIL failed out of $TOTAL"
echo "  Finished: $(date)"
echo "  Receipt log: $RECEIPT_LOG"
echo ""
echo "Next steps:"
echo "  1. Check receipt log for any failures"
echo "  2. Run 05_collect_accessions.py to gather GCA accessions"
echo "  3. Accessions may take 24-48h to appear in ENA browser"
echo "============================================================"
