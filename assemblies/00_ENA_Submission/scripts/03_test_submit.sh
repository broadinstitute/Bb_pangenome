#!/usr/bin/env bash
#
# 03_test_submit.sh — Submit a few assemblies to ENA's TEST server.
#
# Usage:
#   ./03_test_submit.sh <manifest_dir> <webin_user> <webin_password> [n_test] [webin_cli_jar]
#
# Submits the first N manifests (default 3) to ENA's test environment.
# This does NOT create real accessions.

set -euo pipefail

MANIFEST_DIR="${1:?Usage: $0 <manifest_dir> <webin_user> [n_test] [webin_cli_jar]}"
WEBIN_USER="${2:?Provide Webin username}"
#why is this plaintext???
#WEBIN_PASS="${3:?Provide Webin password}"
stty -echo
printf "Webin password: "
read WEBIN_PASS
stty echo
echo
N_TEST="${3:-3}"
WEBIN_CLI="${4:-webin-cli.jar}"
LOG_DIR="${MANIFEST_DIR}/logs/test_submit"

mkdir -p "$LOG_DIR"

COUNT=0
PASS=0
FAIL=0

echo "============================================================"
echo "TEST submission: first $N_TEST assemblies"
echo "  Server: ENA TEST (no real accessions created)"
echo "============================================================"
echo ""

for manifest in "$MANIFEST_DIR"/manifest_*.txt; do
    [ -f "$manifest" ] || continue
    [ "$COUNT" -ge "$N_TEST" ] && break
    COUNT=$((COUNT + 1))

    sample=$(basename "$manifest" .txt | sed 's/^manifest_//')
    logfile="$LOG_DIR/${sample}.log"

    echo -n "  [$COUNT/$N_TEST] Submitting: $sample ... "
    mkdir -p "$LOG_DIR/$sample"
    if java -jar "$WEBIN_CLI" \
        -context genome \
        -userName "$WEBIN_USER" \
        -password "$WEBIN_PASS" \
        -manifest "$manifest" \
        -outputDir "$LOG_DIR/$sample/" \
        -test \
        -submit \
        > "$logfile" 2>&1; then
        echo "OK"
        PASS=$((PASS + 1))
    else
        echo "FAIL (see $logfile)"
        FAIL=$((FAIL + 1))
    fi

    sleep 2
done

echo ""
echo "============================================================"
echo "Test submission: $PASS passed, $FAIL failed out of $COUNT"
echo "Verify at: https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/"
echo "============================================================"
