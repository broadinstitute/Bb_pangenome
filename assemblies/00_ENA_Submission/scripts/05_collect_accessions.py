#!/usr/bin/env python3
"""
Collect assembly accessions from ENA submission logs.

Parses webin-cli output logs to extract ERZ receipt accessions,
then queries ENA API to get GCA accessions when available.

Usage:
    python 05_collect_accessions.py \
        --log-dir manifests/logs/submit/ \
        --receipt-log manifests/submission_receipts_*.tsv \
        --output accessions.tsv
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.request
from pathlib import Path


def parse_webin_log(log_path: Path) -> dict:
    """Extract accession info from a webin-cli submission log."""
    result = {
        "log_file": str(log_path),
        "erz_accession": None,
        "status": "unknown",
        "errors": [],
    }

    text = log_path.read_text(errors="replace")

    # Look for ERZ accession in the log
    erz_match = re.search(r"(ERZ\d+)", text)
    if erz_match:
        result["erz_accession"] = erz_match.group(1)
        result["status"] = "submitted"

    # Look for success indicators
    if "has been completed" in text.lower() or "successfully" in text.lower():
        result["status"] = "submitted"

    # Look for errors
    if "error" in text.lower():
        error_lines = [line.strip() for line in text.splitlines()
                      if "error" in line.lower() and line.strip()]
        result["errors"] = error_lines[:5]
        if not result["erz_accession"]:
            result["status"] = "failed"

    return result


def query_ena_accession(erz_accession: str) -> str | None:
    """Query ENA API to get the GCA accession for an ERZ receipt."""
    url = f"https://www.ebi.ac.uk/ena/submit/report/analyses/{erz_accession}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            if data and isinstance(data, list) and len(data) > 0:
                return data[0].get("accession")
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser(description="Collect ENA assembly accessions.")
    ap.add_argument("--log-dir", type=Path, required=True,
                    help="Directory containing webin-cli submission logs")
    ap.add_argument("--receipt-log", type=Path, default=None,
                    help="Receipt TSV from batch submission script")
    ap.add_argument("--output", type=Path, default=Path("accessions.tsv"),
                    help="Output accession mapping file")
    ap.add_argument("--query-ena", action="store_true",
                    help="Query ENA API for GCA accessions (may not be ready yet)")
    args = ap.parse_args()

    results = []

    # Parse all log files
    log_files = sorted(args.log_dir.rglob("*.log"))
    print(f"Found {len(log_files)} log files in {args.log_dir}")

    for log_path in log_files:
        # Extract assembly name from log filename or parent dir
        sample = log_path.stem
        if sample == "submit" or sample == "validate":
            continue

        info = parse_webin_log(log_path)
        info["assembly_name"] = sample
        info["gca_accession"] = None

        # Optionally query ENA for GCA
        if args.query_ena and info["erz_accession"]:
            print(f"  Querying ENA for {info['erz_accession']}...", end=" ")
            gca = query_ena_accession(info["erz_accession"])
            if gca:
                info["gca_accession"] = gca
                print(f"-> {gca}")
            else:
                print("not yet available")
            time.sleep(1)  # rate limit

        results.append(info)

    # If we have a receipt log, merge with it
    if args.receipt_log and args.receipt_log.exists():
        receipt_data = {}
        with open(args.receipt_log) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                receipt_data[row.get("assembly_name", "")] = row

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "assembly_name", "erz_accession", "gca_accession",
            "status", "errors"
        ])
        for r in results:
            writer.writerow([
                r["assembly_name"],
                r.get("erz_accession", ""),
                r.get("gca_accession", ""),
                r["status"],
                "; ".join(r.get("errors", []))
            ])

    print(f"\nAccession mapping written to {args.output}")

    # Summary
    submitted = sum(1 for r in results if r["status"] == "submitted")
    failed = sum(1 for r in results if r["status"] == "failed")
    with_gca = sum(1 for r in results if r.get("gca_accession"))
    print(f"  Submitted: {submitted}")
    print(f"  Failed: {failed}")
    if args.query_ena:
        print(f"  GCA assigned: {with_gca}")
        if with_gca < submitted:
            print(f"  (remaining may take 24-48h to appear)")


if __name__ == "__main__":
    main()
