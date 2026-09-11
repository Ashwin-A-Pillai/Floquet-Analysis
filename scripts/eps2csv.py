#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse a single EPS spectrum file and export two CSVs:
  - csv_imag/eps_imag.csv  (E, Im(eps))
  - csv_real/eps_real.csv  (E, Re(eps))

Usage:
  python eps2csv.py -I o-IPA.eps_q1_ip
  python eps2csv.py -I o-IPA.eps_q1_ip --out-prefix o-IPA

Python 3.6+, no third-party deps.
"""

from __future__ import print_function
import os
import sys
import csv
import argparse


def parse_eps_file(path):
    """
    Returns: list of (E, im_eps, re_eps) as floats.
    Skips header lines starting with '#'.
    Expects at least 3 columns:
      E[eV], Im(eps), Re(eps), ...
    """
    data = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                e = float(parts[0])
                imv = float(parts[1])
                rev = float(parts[2])
            except ValueError:
                continue
            data.append((e, imv, rev))
    return data


def ensure_dir(d):
    if d and not os.path.isdir(d):
        os.makedirs(d)


def main(argv):
    p = argparse.ArgumentParser(description="Convert EPS text file to imag/real CSVs.")
    p.add_argument("-I", "--input", required=True, help="Input EPS file (e.g. o-IPA.eps_q1_ip)")
    p.add_argument("--out-imag-dir", default="csv_imag", help="Output directory for imag CSV")
    p.add_argument("--out-real-dir", default="csv_real", help="Output directory for real CSV")
    p.add_argument("--out-prefix", default="eps", help="Output filename prefix (default: eps)")
    args = p.parse_args(argv[1:])

    infile = args.input
    if not os.path.isfile(infile):
        print("ERROR: input file not found: {}".format(infile), file=sys.stderr)
        return 1

    rows = parse_eps_file(infile)
    if not rows:
        print("ERROR: no data parsed from: {}".format(infile), file=sys.stderr)
        return 1

    # Sort by energy, stable
    rows.sort(key=lambda t: t[0])

    ensure_dir(args.out_imag_dir)
    ensure_dir(args.out_real_dir)

    out_imag_path = os.path.join(args.out_imag_dir, "{}_imag.csv".format(args.out_prefix))
    out_real_path = os.path.join(args.out_real_dir, "{}_real.csv".format(args.out_prefix))

    # Write IMAG
    with open(out_imag_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E[1] [eV]", "Im(eps)"])
        for e, imv, _ in rows:
            w.writerow(["{:.6f}".format(e), "{:.8f}".format(imv)])

    # Write REAL
    with open(out_real_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E[1] [eV]", "Re(eps)"])
        for e, _, rev in rows:
            w.writerow(["{:.6f}".format(e), "{:.8f}".format(rev)])

    print("Wrote:\n  {}\n  {}".format(out_imag_path, out_real_path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))