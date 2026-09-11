#!/usr/bin/env python3
# compute_rms_neighbors.py  (Python 3.6.8)

import csv
import math
import re
import sys
from typing import List, Tuple

# Match either "Re(eps)_###" or "Im(eps)_###"
EPS_COL_RE = re.compile(r'^(Re|Im)\(eps\)_(\d+)$')

def find_eps_cols(header: List[str]) -> Tuple[List[Tuple[int, str]], str]:
    """
    Return (cols, series_kind)
      - cols: list of (index, label) for columns named like 'Re(eps)_###' or 'Im(eps)_###',
              preserving their order in the CSV.
      - series_kind: 'Re' or 'Im' (all epsilon columns must be the same kind).
    """
    cols = []
    series_kind = None
    for idx, name in enumerate(header):
        m = EPS_COL_RE.match(name.strip())
        if m:
            kind, label = m.group(1), m.group(2)  # kind: 'Re' or 'Im'; label: '300', '4000', ...
            if series_kind is None:
                series_kind = kind
            elif kind != series_kind:
                raise ValueError("Mixed 'Re(eps)_###' and 'Im(eps)_###' columns in one file; "
                                 "please use a single kind per CSV.")
            cols.append((idx, label))
    if len(cols) < 2:
        raise ValueError("Need at least two epsilon columns named like 'Re(eps)_###' or 'Im(eps)_###'.")
    return cols, series_kind

def compute_rms_neighbors(input_path: str, output_path: str) -> None:
    with open(input_path, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)

        eps_cols, series_kind = find_eps_cols(header)  # e.g., [(1,'200'), (2,'300'), ...], 'Re' or 'Im'
        num_pairs = len(eps_cols) - 1
        if num_pairs < 1:
            raise ValueError("Not enough neighboring columns to compare.")

        # Accumulate sum of squared diffs per neighbor pair
        ssd = [0.0] * num_pairs
        N = 0

        for row in reader:
            if not row or all(c.strip() == "" for c in row):
                continue  # skip blank lines

            # Convert only the epsilon columns for this row to floats
            vals = []
            for idx, _ in eps_cols:
                try:
                    vals.append(float(row[idx]))
                except (ValueError, IndexError):
                    # If any value is missing/invalid, skip this row entirely
                    vals = None
                    break
            if vals is None:
                continue

            # Add squared differences for each neighboring pair
            for p in range(num_pairs):
                d = vals[p] - vals[p + 1]
                ssd[p] += d * d
            N += 1

        if N == 0:
            raise ValueError("No valid data rows were found.")

        # Prepare results: label is the *second* column in each pair
        results = []
        for p in range(num_pairs):
            label = eps_cols[p + 1][1]  # e.g., '300' for (200 vs 300) or '4000' for (3000 vs 4000)
            rms = math.sqrt(ssd[p] / N)
            results.append((label, rms))

    # Write output CSV
    with open(output_path, 'w', newline='', encoding='utf-8') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(['label', 'rms_between_neighbors'])  # label is numeric string like '300' or '4000'
        for label, rms in results:
            writer.writerow([label, "{:.9f}".format(rms)])

if __name__ == '__main__':
    # Usage: python compute_rms_neighbors.py input.csv [output.csv]
    if len(sys.argv) < 2:
        print("Usage: python compute_rms_neighbors.py <input.csv> [output.csv]")
        sys.exit(1)
    input_csv = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else 'rms_between_neighbors.csv'
    compute_rms_neighbors(input_csv, output_csv)
    print("Wrote:", output_csv)
