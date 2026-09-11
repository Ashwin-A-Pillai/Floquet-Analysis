#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eps_rms.py

Generic EPS RMS convergence script.

This replaces both:

    eps_rms.py
    eps_rms_nband.py

It keeps the eps_rms_nband.py command-line interface:

    --nbands
    --pattern
    --out-dir
    --no-plot

The name --nbands is retained for compatibility, but the values can represent
any convergence variable, for example:

    nbnd values: 70 80 90 100
    nk values  : 6 8 10 12 14 16 18

Input pattern must use {nbnd} as placeholder, for example:

    eps[{nbnd}].out
    o-IPA-{nbnd}.eps_q1_inv_rpa_dyson

Additional option:

    --x-label

Examples:

    # nbnd convergence
    python3 eps_rms.py \
      --nbands 90 100 110 120 140 160 \
      --pattern "eps[{nbnd}].out" \
      --out-dir eps_nbnd/rms \
      --x-label "Number of bands, nbnd"

    # nk convergence
    python3 eps_rms.py \
      --nbands 6 8 10 12 14 16 18 \
      --pattern "eps[{nbnd}].out" \
      --out-dir eps_nk/rms \
      --x-label "nk"

Reads EPS files with columns:

    E[eV], Im(eps), Re(eps), ...

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


def get_common_energies(data_by_param, params):
    common = None

    for param in params:
        keys = set(data_by_param[param].keys())

        if common is None:
            common = keys
        else:
            common = common.intersection(keys)

    if not common:
        raise ValueError("No common energy values found across all EPS files.")

    return sorted(common, key=lambda x: float(x))


def write_combined_csv(output_path, data_by_param, params, energies, kind, column_label):
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

        for param in params:
            header.append("{}_{}_{}".format(col_base, column_label, param))

        writer.writerow(header)

        for e_key in energies:
            row = [e_key]

            for param in params:
                value = data_by_param[param][e_key][value_index]
                row.append("{:.8f}".format(value))

            writer.writerow(row)


def get_values(data_by_param, param, energies, kind):
    if kind == "real":
        value_index = 2
    elif kind == "imag":
        value_index = 1
    else:
        raise ValueError("kind must be 'real' or 'imag'")

    vals = []

    for e_key in energies:
        vals.append(data_by_param[param][e_key][value_index])

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


def compute_neighbor_rms(data_by_param, params, energies, kind):
    """
    Returns:
        [(param_1, param_2, rms), ...]
    """
    results = []

    for param1, param2 in zip(params[:-1], params[1:]):
        vals1 = get_values(data_by_param, param1, energies, kind)
        vals2 = get_values(data_by_param, param2, energies, kind)

        r = rms_between_arrays(vals1, vals2)
        results.append((param1, param2, r))

    return results


def write_neighbor_rms(output_path, neighbor_results, n_energy, column_label):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "{}_1".format(column_label),
            "{}_2".format(column_label),
            "plot_x_{}".format(column_label),
            "N_common_energy_points",
            "rms",
        ])

        for param1, param2, r in neighbor_results:
            writer.writerow([
                param1,
                param2,
                param2,
                n_energy,
                "{:.9f}".format(r),
            ])


def write_all_pair_rms(output_path, data_by_param, params, energies, kind, column_label):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "{}_1".format(column_label),
            "{}_2".format(column_label),
            "N_common_energy_points",
            "rms",
        ])

        for param1, param2 in itertools.combinations(params, 2):
            vals1 = get_values(data_by_param, param1, energies, kind)
            vals2 = get_values(data_by_param, param2, energies, kind)

            r = rms_between_arrays(vals1, vals2)

            writer.writerow([
                param1,
                param2,
                len(energies),
                "{:.9f}".format(r),
            ])


def plot_neighbor_rms(output_path, neighbor_results, title, ylabel, x_label):
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

    for param1, param2, r in neighbor_results:
        x.append(float(param2))
        y.append(r)
        labels.append("({},{})".format(param1, param2))

    plt.figure()
    plt.plot(x, y, marker="o")
    plt.xlabel(x_label)
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


def plot_real_imag_together(output_path, real_results, imag_results, x_label):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not found. Skipping plot:", output_path)
        return

    x_real = []
    y_real = []

    for param1, param2, r in real_results:
        x_real.append(float(param2))
        y_real.append(r)

    x_imag = []
    y_imag = []

    for param1, param2, r in imag_results:
        x_imag.append(float(param2))
        y_imag.append(r)

    plt.figure()
    plt.plot(x_real, y_real, marker="o", label="Re(eps)")
    plt.plot(x_imag, y_imag, marker="s", label="Im(eps)")
    plt.xlabel(x_label)
    plt.ylabel("RMS between neighboring values")
    plt.title("Neighbor RMS convergence of Re(eps) and Im(eps)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main(argv):
    parser = argparse.ArgumentParser(
        description="Parse EPS files for convergence values, compute RMS, and plot convergence."
    )

    parser.add_argument(
        "--nbands",
        nargs="+",
        default=["70", "80", "90", "100", "110", "120"],
        help=(
            "Convergence values. The name is kept for compatibility with "
            "eps_rms_nband.py. These may be nbnd, nk, or any other parameter. "
            "Default: 70 80 90 100 110 120"
        ),
    )

    parser.add_argument(
        "--pattern",
        default="eps[{nbnd}].out",
        help=(
            "Input EPS filename pattern. Use {nbnd} as placeholder. "
            "Default: eps[{nbnd}].out"
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory. Default: current directory.",
    )

    parser.add_argument(
        "--x-label",
        default="Number of bands, nbnd",
        help=(
            "Label for the x-axis in plots. "
            "Examples: 'Number of bands, nbnd', 'nk'. "
            "Default: 'Number of bands, nbnd'"
        ),
    )

    parser.add_argument(
        "--column-label",
        default="nbnd",
        help=(
            "Label used in CSV column names. "
            "Examples: nbnd, nk. Default: nbnd"
        ),
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Compute CSV files only; do not generate PNG plots.",
    )

    args = parser.parse_args(argv[1:])

    params = args.nbands
    ensure_dir(args.out_dir)

    data_by_param = {}

    print("Reading EPS files:")

    for param in params:
        path = args.pattern.format(nbnd=param)

        if not os.path.isfile(path):
            print("ERROR: file not found: {}".format(path), file=sys.stderr)
            return 1

        data = parse_eps_file(path)

        if not data:
            print("ERROR: no data parsed from: {}".format(path), file=sys.stderr)
            return 1

        data_by_param[param] = data

        print(
            "  {} = {:>4s} : {} points : {}".format(
                args.column_label,
                param,
                len(data),
                path,
            )
        )

    energies = get_common_energies(data_by_param, params)

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

    write_combined_csv(
        combined_real,
        data_by_param,
        params,
        energies,
        "real",
        args.column_label,
    )

    write_combined_csv(
        combined_imag,
        data_by_param,
        params,
        energies,
        "imag",
        args.column_label,
    )

    real_neighbor_results = compute_neighbor_rms(
        data_by_param,
        params,
        energies,
        "real",
    )

    imag_neighbor_results = compute_neighbor_rms(
        data_by_param,
        params,
        energies,
        "imag",
    )

    write_neighbor_rms(
        rms_real_neighbors,
        real_neighbor_results,
        len(energies),
        args.column_label,
    )

    write_neighbor_rms(
        rms_imag_neighbors,
        imag_neighbor_results,
        len(energies),
        args.column_label,
    )

    write_all_pair_rms(
        rms_real_all_pairs,
        data_by_param,
        params,
        energies,
        "real",
        args.column_label,
    )

    write_all_pair_rms(
        rms_imag_all_pairs,
        data_by_param,
        params,
        energies,
        "imag",
        args.column_label,
    )

    if not args.no_plot:
        plot_neighbor_rms(
            plot_real,
            real_neighbor_results,
            "Neighbor RMS convergence of Re(eps)",
            "RMS of Re(eps)",
            args.x_label,
        )

        plot_neighbor_rms(
            plot_imag,
            imag_neighbor_results,
            "Neighbor RMS convergence of Im(eps)",
            "RMS of Im(eps)",
            args.x_label,
        )

        plot_real_imag_together(
            plot_both,
            real_neighbor_results,
            imag_neighbor_results,
            args.x_label,
        )

    print("")
    print("Neighbor RMS values for Re(eps):")
    for param1, param2, r in real_neighbor_results:
        print("  x = {:>4s}, RMS({},{}) = {:.9f}".format(param2, param1, param2, r))

    print("")
    print("Neighbor RMS values for Im(eps):")
    for param1, param2, r in imag_neighbor_results:
        print("  x = {:>4s}, RMS({},{}) = {:.9f}".format(param2, param1, param2, r))

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