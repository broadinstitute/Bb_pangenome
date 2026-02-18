#!/usr/bin/env bash
#
# 02_validate.sh — Run webin-cli validation on all generated manifests.
#
# Usage:
#   ./02_validate.sh <manifest_dir> <webin_user> <webin_password> [webin_cli_jar]
#
# Validates each manifest without submitting. Reports pass/fail for each.

set -euo pipefail

MANIFEST_DIR="${1:?Usage: $0 <manifest_dir> <webin_user> [webin_cli_jar]}"
WEBIN_USER="${2:?Provide Webin username}"
#why is this plaintext???
#WEBIN_PASS="${3:?Provide Webin password}"
stty -echo
printf "Webin password: "
read WEBIN_PASS
stty echo
echo
WEBIN_CLI="${3:-webin-cli.jar}"
LOG_DIR="${MANIFEST_DIR}/logs/validate"
mkdir -p "$LOG_DIR"

# Check webin-cli exists
if [ ! -f "$WEBIN_CLI" ]; then
    echo "ERROR: webin-cli.jar not found at $WEBIN_CLI"
    echo "Download from: https://github.com/enaDataSubmission/webin-cli/releases"
    exit 1
fi

# Check java
if ! command -v java &>/dev/null; then
    echo "ERROR: java not found. webin-cli requires Java 11+."
    exit 1
fi

PASS=0
FAIL=0
TOTAL=0

echo "============================================================"
echo "Validating manifests in: $MANIFEST_DIR"
echo "============================================================"
echo ""

for manifest in "$MANIFEST_DIR"/manifest_*.txt; do
    [ -f "$manifest" ] || continue
    TOTAL=$((TOTAL + 1))

    sample=$(basename "$manifest" .txt | sed 's/^manifest_//')
    logfile="$LOG_DIR/${sample}.log"

    echo -n "  Validating: $sample ... "
    mkdir -p "$LOG_DIR/$sample"
    if java -jar "$WEBIN_CLI" \
        -context genome \
        -userName "$WEBIN_USER" \
        -password "$WEBIN_PASS" \
        -manifest "$manifest" \
        -outputDir "$(realpath "$LOG_DIR/$sample/")" \
        -validate \
        > "$logfile" 2>&1; then
        echo "PASS"
        PASS=$((PASS + 1))
    else
        echo "FAIL (see $logfile)"
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "============================================================"
echo "Validation complete: $PASS passed, $FAIL failed, $TOTAL total"
echo "Logs: $LOG_DIR"
echo "============================================================"

if [ "$FAIL" -gt 0 ]; then
    echo ""
    echo "Failed manifests:"
    for manifest in "$MANIFEST_DIR"/manifest_*.txt; do
        [ -f "$manifest" ] || continue
        sample=$(basename "$manifest" .txt | sed 's/^manifest_//')
        logfile="$LOG_DIR/${sample}.log"
        if grep -qi "error" "$logfile" 2>/dev/null; then
            echo "  $sample:"
            grep -i "error" "$logfile" | head -5 | sed 's/^/    /'
        fi
    done
    exit 1
fi
