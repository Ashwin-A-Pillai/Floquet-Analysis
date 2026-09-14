#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eps_rms_nbnd.py

Reads EPS files for different nbnd values, for example:

    eps[70].out
    eps[80].out
    eps[90].out
    eps[100].out
    eps[110].out
    eps[120].out

Extracts:
    E[eV], Im(eps), Re(eps)

Writes:
    combined_real.csv
    combined_imag.csv

    rms_real_neighbors.csv
    rms_imag_neighbors.csv

    rms_real_all_pairs.csv
    rms_imag_all_pairs.csv

Plots:
    rms_real_neighbors.png
    rms_imag_neighbors.png
    rms_real_imag_neighbors.png

For neighboring RMS:
    RMS(70,80)   is plotted at x = 80
    RMS(80,90)   is plotted at x = 90
    RMS(90,100)  is plotted at x = 100
"""

from __future__ import print_function

import os
import sys
import csv
import math
import argparse
import itertools


def parse_eps_file(path):
    """
    Expected columns:
        E[eV], Im(eps), Re(eps), ...

    Returns:
        energy_key -> (energy_float, im_eps, re_eps)
    """
    data = {}

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

            e_key = "{:.6f}".format(e)
            data[e_key] = (e, imv, rev)

    return data


def ensure_dir(path):
    if path and path != "." and not os.path.isdir(path):
        os.makedirs(path)


def get_common_energies(data_by_nbnd, nbands):
    common = None

    for nbnd in nbands:
        keys = set(data_by_nbnd[nbnd].keys())

        if common is None:
            common = keys
        else:
            common = common.intersection(keys)

    if not common:
        raise ValueError("No common energy values found across all nbnd files.")

    return sorted(common, key=lambda x: float(x))


def write_combined_csv(output_path, data_by_nbnd, nbands, energies, kind):
    if kind == "real":
        col_base = "Re(eps)"
        value_index = 2
    elif kind == "imag":
        col_base = "Im(eps)"
        value_index = 1
    else:
        raise ValueError("kind must be 'real' or 'imag'")

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        header = ["E[1] [eV]"]

        for nbnd in nbands:
            header.append("{}_nbnd_{}".format(col_base, nbnd))

        writer.writerow(header)

        for e_key in energies:
            row = [e_key]

            for nbnd in nbands:
                value = data_by_nbnd[nbnd][e_key][value_index]
                row.append("{:.8f}".format(value))

            writer.writerow(row)


def get_values(data_by_nbnd, nbnd, energies, kind):
    if kind == "real":
        value_index = 2
    elif kind == "imag":
        value_index = 1
    else:
        raise ValueError("kind must be 'real' or 'imag'")

    vals = []

    for e_key in energies:
        vals.append(data_by_nbnd[nbnd][e_key][value_index])

    return vals


def rms_between_arrays(vals1, vals2):
    if len(vals1) != len(vals2):
        raise ValueError("Cannot compute RMS for arrays of different length.")

    if len(vals1) == 0:
        raise ValueError("Cannot compute RMS for empty arrays.")

    ssd = 0.0

    for a, b in zip(vals1, vals2):
        d = a - b
        ssd += d * d

    return math.sqrt(ssd / len(vals1))


def compute_neighbor_rms(data_by_nbnd, nbands, energies, kind):
    """
    Returns:
        [(nbnd_1, nbnd_2, rms), ...]
    """
    results = []

    for nbnd1, nbnd2 in zip(nbands[:-1], nbands[1:]):
        vals1 = get_values(data_by_nbnd, nbnd1, energies, kind)
        vals2 = get_values(data_by_nbnd, nbnd2, energies, kind)

        r = rms_between_arrays(vals1, vals2)
        results.append((nbnd1, nbnd2, r))

    return results


def write_neighbor_rms(output_path, neighbor_results, n_energy):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "nbnd_1",
            "nbnd_2",
            "plot_x_nbnd",
            "N_common_energy_points",
            "rms",
        ])

        for nbnd1, nbnd2, r in neighbor_results:
            writer.writerow([
                nbnd1,
                nbnd2,
                nbnd2,
                n_energy,
                "{:.9f}".format(r),
            ])


def write_all_pair_rms(output_path, data_by_nbnd, nbands, energies, kind):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "nbnd_1",
            "nbnd_2",
            "N_common_energy_points",
            "rms",
        ])

        for nbnd1, nbnd2 in itertools.combinations(nbands, 2):
            vals1 = get_values(data_by_nbnd, nbnd1, energies, kind)
            vals2 = get_values(data_by_nbnd, nbnd2, energies, kind)

            r = rms_between_arrays(vals1, vals2)

            writer.writerow([
                nbnd1,
                nbnd2,
                len(energies),
                "{:.9f}".format(r),
            ])


def plot_neighbor_rms(output_path, neighbor_results, title, ylabel):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not found. Skipping plot:", output_path)
        return

    x = []
    y = []
    labels = []

    for nbnd1, nbnd2, r in neighbor_results:
        x.append(float(nbnd2))
        y.append(r)
        labels.append("({},{})".format(nbnd1, nbnd2))

    plt.figure()
    plt.plot(x, y, marker="o")
    plt.xlabel("Number of bands, nbnd")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)

    for xi, yi, label in zip(x, y, labels):
        plt.annotate(
            label,
            (xi, yi),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_real_imag_together(output_path, real_results, imag_results):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not found. Skipping plot:", output_path)
        return

    x_real = []
    y_real = []

    for nbnd1, nbnd2, r in real_results:
        x_real.append(float(nbnd2))
        y_real.append(r)

    x_imag = []
    y_imag = []

    for nbnd1, nbnd2, r in imag_results:
        x_imag.append(float(nbnd2))
        y_imag.append(r)

    plt.figure()
    plt.plot(x_real, y_real, marker="o", label="Re(eps)")
    plt.plot(x_imag, y_imag, marker="s", label="Im(eps)")
    plt.xlabel("Number of bands, nbnd")
    plt.ylabel("RMS between neighboring nbnd values")
    plt.title("Neighbor RMS convergence of Re(eps) and Im(eps) with nbnd")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main(argv):
    parser = argparse.ArgumentParser(
        description="Parse EPS files for several nbnd values, compute RMS, and plot convergence."
    )

    parser.add_argument(
        "--nbands",
        nargs="+",
        default=["70", "80", "90", "100", "110", "120"],
        help="nbnd values. Default: 70 80 90 100 110 120",
    )

    parser.add_argument(
        "--pattern",
        default="eps[{nbnd}].out",
        help="Input EPS filename pattern. Use {nbnd} as placeholder. Default: eps[{nbnd}].out",
    )

    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory. Default: current directory.",
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Compute CSV files only; do not generate PNG plots.",
    )

    args = parser.parse_args(argv[1:])

    nbands = args.nbands
    ensure_dir(args.out_dir)

    data_by_nbnd = {}

    print("Reading EPS files:")

    for nbnd in nbands:
        path = args.pattern.format(nbnd=nbnd)

        if not os.path.isfile(path):
            print("ERROR: file not found: {}".format(path), file=sys.stderr)
            return 1

        data = parse_eps_file(path)

        if not data:
            print("ERROR: no data parsed from: {}".format(path), file=sys.stderr)
            return 1

        data_by_nbnd[nbnd] = data

        print("  nbnd = {:>4s} : {} points : {}".format(nbnd, len(data), path))

    energies = get_common_energies(data_by_nbnd, nbands)

    print("")
    print("Using {} common energy points for RMS.".format(len(energies)))

    combined_real = os.path.join(args.out_dir, "combined_real.csv")
    combined_imag = os.path.join(args.out_dir, "combined_imag.csv")

    rms_real_neighbors = os.path.join(args.out_dir, "rms_real_neighbors.csv")
    rms_imag_neighbors = os.path.join(args.out_dir, "rms_imag_neighbors.csv")

    rms_real_all_pairs = os.path.join(args.out_dir, "rms_real_all_pairs.csv")
    rms_imag_all_pairs = os.path.join(args.out_dir, "rms_imag_all_pairs.csv")

    plot_real = os.path.join(args.out_dir, "rms_real_neighbors.png")
    plot_imag = os.path.join(args.out_dir, "rms_imag_neighbors.png")
    plot_both = os.path.join(args.out_dir, "rms_real_imag_neighbors.png")

    write_combined_csv(combined_real, data_by_nbnd, nbands, energies, "real")
    write_combined_csv(combined_imag, data_by_nbnd, nbands, energies, "imag")

    real_neighbor_results = compute_neighbor_rms(data_by_nbnd, nbands, energies, "real")
    imag_neighbor_results = compute_neighbor_rms(data_by_nbnd, nbands, energies, "imag")

    write_neighbor_rms(rms_real_neighbors, real_neighbor_results, len(energies))
    write_neighbor_rms(rms_imag_neighbors, imag_neighbor_results, len(energies))

    write_all_pair_rms(rms_real_all_pairs, data_by_nbnd, nbands, energies, "real")
    write_all_pair_rms(rms_imag_all_pairs, data_by_nbnd, nbands, energies, "imag")

    if not args.no_plot:
        plot_neighbor_rms(
            plot_real,
            real_neighbor_results,
            "Neighbor RMS convergence of Re(eps) with nbnd",
            "RMS of Re(eps)",
        )

        plot_neighbor_rms(
            plot_imag,
            imag_neighbor_results,
            "Neighbor RMS convergence of Im(eps) with nbnd",
            "RMS of Im(eps)",
        )

        plot_real_imag_together(
            plot_both,
            real_neighbor_results,
            imag_neighbor_results,
        )

    print("")
    print("Neighbor RMS values for Re(eps):")
    for nbnd1, nbnd2, r in real_neighbor_results:
        print("  x = {:>4s}, RMS({},{}) = {:.9f}".format(nbnd2, nbnd1, nbnd2, r))

    print("")
    print("Neighbor RMS values for Im(eps):")
    for nbnd1, nbnd2, r in imag_neighbor_results:
        print("  x = {:>4s}, RMS({},{}) = {:.9f}".format(nbnd2, nbnd1, nbnd2, r))

    print("")
    print("Wrote:")
    print("  {}".format(combined_real))
    print("  {}".format(combined_imag))
    print("  {}".format(rms_real_neighbors))
    print("  {}".format(rms_imag_neighbors))
    print("  {}".format(rms_real_all_pairs))
    print("  {}".format(rms_imag_all_pairs))

    if not args.no_plot:
        print("  {}".format(plot_real))
        print("  {}".format(plot_imag))
        print("  {}".format(plot_both))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))