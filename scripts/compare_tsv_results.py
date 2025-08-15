#!/usr/bin/env python3
import sys, csv
from pathlib import Path
from collections import Counter

def read_tsv_by_key(path):
    d = {}
    seen = Counter()
    with Path(path).open(newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if not row:
                continue
            key, *vals = row
            seen[key] += 1
            d[key] = vals
    dups = [k for k,c in seen.items() if c>1]
    if dups:
        sys.stderr.write(f"warning: {path} has duplicate keys: {', '.join(dups[:10])}"
                         f"{' ...' if len(dups)>10 else ''}\n")
    return d

def main(f1, f2, identical_file, different_file):
    a = read_tsv_by_key(f1)
    b = read_tsv_by_key(f2)

    keys_a = set(a)
    keys_b = set(b)
    both = keys_a & keys_b

    identical = []
    differences = []

    # Compare overlapping keys
    for k in sorted(both):
        if a[k] == b[k]:
            identical.append(k)
        else:
            differences.append((k, "mismatch"))

    # Keys unique to file1
    for k in sorted(keys_a - keys_b):
        differences.append((k, f"only in {Path(f1).name}"))

    # Keys unique to file2
    for k in sorted(keys_b - keys_a):
        differences.append((k, f"only in {Path(f2).name}"))

    # Write identical.txt
    with open(identical_file, "w") as f:
        for k in identical:
            f.write(f"{k}\n")

    # Write different.txt with source info
    with open(different_file, "w") as f:
        for k, reason in differences:
            f.write(f"{k}\t{reason}\n")

    # Summary
    sys.stderr.write(
        f"{f1}: {len(keys_a)} rows, {f2}: {len(keys_b)} rows\n"
        f"Overlap: {len(both)} | Identical: {len(identical)} | "
        f"Different: {len(differences)}\n"
    )

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write(
            f"usage: {Path(sys.argv[0]).name} file1.tsv file2.tsv identical.txt different.txt\n"
        )
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
