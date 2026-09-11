#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qepy.lattice import Path as QEPath

from yambopy import (
    YamboLatticeDB,
    YamboVbandsDB,
    findallk_qe,
    get_qebands_interpolate,
    plot_2D_kdist,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="All-k Floquet analysis and interpolated band plotting."
    )

    parser.add_argument(
        "--workdir",
        default=".",
        help="Yambo/FixSymm working directory.",
    )

    parser.add_argument(
        "--calc",
        default="fl_sample",
        help="Floquet database directory.",
    )

    parser.add_argument(
        "--save-db",
        default="SAVE/ns.db1",
        help="Path to ns.db1.",
    )

    parser.add_argument(
        "--path",
        required=True,
        help=(
            "Band path: "
            "'G:0,0,0;X:0.5,0,0.5;W:0.5,0.25,0.75;...'"
        ),
    )

    parser.add_argument(
        "--npoints",
        type=int,
        default=10,
        help="Interpolation divisions per path segment.",
    )

    parser.add_argument(
        "--segment-divisions",
        default=None,
        help=(
            "Optional comma-separated segment divisions. "
            "Overrides --npoints."
        ),
    )

    parser.add_argument(
        "--eta",
        type=int,
        default=2,
        help="Physical Floquet harmonic eta.",
    )

    parser.add_argument(
        "--initial-band",
        type=int,
        default=4,
        help="1-based physical initial band.",
    )

    parser.add_argument(
        "--projected-band",
        type=int,
        default=5,
        help="1-based physical projected band in the Floquet basis.",
    )

    parser.add_argument(
        "--report-file",
        default="results_eval.dat",
        help="findallk_qe report file.",
    )

    parser.add_argument(
        "--band-pdf",
        default="QE_band_interp.pdf",
        help="Output interpolated band PDF.",
    )

    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="Optional lower energy limit in eV.",
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Optional upper energy limit in eV.",
    )

    parser.add_argument(
        "--plot-kdist",
        action="store_true",
        help="Also generate plot_2D_kdist.",
    )

    return parser.parse_args()


def split_calc_path(calc):
    calc_path = Path(calc)

    if calc_path.is_absolute():
        return str(calc_path.parent), calc_path.name

    parent = calc_path.parent
    folder = "." if str(parent) in ("", ".") else str(parent)

    return folder, calc_path.name


def resolve_path(workdir, value):
    path = Path(value).expanduser()

    if path.is_absolute():
        return path

    return (workdir / path).resolve()


def parse_kpath(spec):
    points = []

    for raw_item in spec.split(";"):
        item = raw_item.strip()

        if not item:
            continue

        if ":" not in item:
            raise ValueError(
                f"Invalid path item '{item}'. "
                "Expected LABEL:kx,ky,kz."
            )

        label, coords_text = item.split(":", 1)

        coords = [
            float(x.strip())
            for x in coords_text.split(",")
        ]

        if len(coords) != 3:
            raise ValueError(
                f"Path point '{item}' must contain exactly 3 coordinates."
            )

        label = label.strip()

        if not label:
            raise ValueError(
                f"Missing label in path item '{item}'."
            )

        points.append([coords, label])

    if len(points) < 2:
        raise ValueError(
            "The k path must contain at least two points."
        )

    return points


def get_segment_divisions(args, nsegments):
    if args.segment_divisions is None:
        return [args.npoints] * nsegments

    divisions = [
        int(x.strip())
        for x in args.segment_divisions.split(",")
        if x.strip()
    ]

    if len(divisions) != nsegments:
        raise ValueError(
            f"--segment-divisions contains {len(divisions)} values, "
            f"but the path has {nsegments} segments."
        )

    if any(x <= 0 for x in divisions):
        raise ValueError(
            "All segment divisions must be positive integers."
        )

    return divisions


def main():
    args = parse_args()

    workdir = Path(args.workdir).expanduser().resolve()

    if not workdir.is_dir():
        raise RuntimeError(
            f"Working directory does not exist: {workdir}"
        )

    os.chdir(workdir)

    internal_save = Path("SAVE/ns.db1")

    if not internal_save.is_file():
        raise RuntimeError(
            f"Missing {workdir / internal_save}. "
            "The fl-analysis VbPP implementation reads SAVE/ns.db1 "
            "relative to --workdir."
        )

    save_db = resolve_path(
        workdir,
        args.save_db,
    )

    if not save_db.is_file():
        raise RuntimeError(
            f"Could not find lattice database: {save_db}"
        )

    calc_folder, calc_name = split_calc_path(args.calc)

    mydb = YamboVbandsDB(
        folder=calc_folder,
        calc=calc_name,
    )

    print(mydb)

    if mydb.fl_order <= 0:
        raise RuntimeError(
            f"Invalid Floquet order in database: {mydb.fl_order}."
        )

    expected_timesteps = 2 * mydb.fl_order + 4

    if mydb.n_timesteps < expected_timesteps:
        raise RuntimeError(
            f"Insufficient Floquet sampling: "
            f"N timesteps={mydb.n_timesteps}, "
            f"FLOrder={mydb.fl_order}, "
            f"expected at least {expected_timesteps}."
        )

    if not (-mydb.fl_order <= args.eta <= mydb.fl_order):
        raise ValueError(
            f"--eta={args.eta} lies outside "
            f"[-{mydb.fl_order}, +{mydb.fl_order}]."
        )

    if not (1 <= args.initial_band <= mydb.n_vbands):
        raise ValueError(
            f"--initial-band={args.initial_band} "
            f"is outside 1..{mydb.n_vbands}."
        )

    basis_first, basis_last = mydb.basis_index

    if not (basis_first <= args.projected_band <= basis_last):
        raise ValueError(
            f"--projected-band={args.projected_band} "
            f"is outside the Floquet basis "
            f"{basis_first}..{basis_last}."
        )

    path_points = parse_kpath(
        args.path
    )

    divisions = get_segment_divisions(
        args,
        len(path_points) - 1,
    )

    print("[INFO] Band path:")
    for coords, label in path_points:
        print(f"  {label:>4s} : {coords}")

    print(
        f"[INFO] Segment divisions: {divisions}"
    )

    qe_path = QEPath(
        path_points,
        divisions,
    )

    fl_eigenvectors, fl_quasienergy, ks_eigenvalues = findallk_qe(
        mydb,
        report_file=args.report_file,
    )

    lat = YamboLatticeDB.from_db_file(
        filename=str(save_db),
        Expand=False,
    )

    ks_bs, fl_bs = get_qebands_interpolate(
        lat,
        qe_path,
        fl_quasienergy,
        ks_eigenvalues,
    )

    fig = plt.figure(
        figsize=(5, 6)
    )

    ax = fig.add_axes(
        [0.16, 0.14, 0.80, 0.80]
    )

    ks_kwargs = dict(
        legend=True,
        c_bands="r",
        label="KS",
    )

    fl_kwargs = dict(
        legend=True,
        c_bands="b",
        linestyle="dashed",
        label="FL",
    )

    if args.ymin is not None and args.ymax is not None:
        ks_kwargs["ylim"] = (
            args.ymin,
            args.ymax,
        )
        fl_kwargs["ylim"] = (
            args.ymin,
            args.ymax,
        )

    ks_bs.plot_ax(
        ax,
        **ks_kwargs,
    )

    fl_bs.plot_ax(
        ax,
        **fl_kwargs,
    )

    plt.savefig(
        args.band_pdf,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"[OK] Wrote interpolated band plot: "
        f"{args.band_pdf}"
    )

    # eta axis is ordered -M,...,+M
    floq_mode = args.eta + mydb.fl_order

    # initial physical band -> zero-based second axis
    floq_initial_index = args.initial_band - 1

    # projected physical band -> zero-based Floquet-basis column
    floq_projected_index = (
        args.projected_band - basis_first
    )

    print(
        f"[INFO] eta={args.eta} "
        f"-> Floquet-mode index {floq_mode}"
    )

    print(
        f"[INFO] initial band {args.initial_band} "
        f"-> array index {floq_initial_index}"
    )

    print(
        f"[INFO] projected band {args.projected_band} "
        f"-> Floquet-basis column "
        f"{floq_projected_index}"
    )

    data = np.abs(
        fl_eigenvectors[
            :,
            floq_initial_index,
            floq_mode,
            floq_projected_index,
        ]
    )

    coeff_file = (
        f"floquet_coeff_eta_{args.eta:+d}_"
        f"b{args.initial_band}_to_b"
        f"{args.projected_band}.dat"
    )

    np.savetxt(
        coeff_file,
        np.column_stack(
            [
                np.arange(
                    1,
                    mydb.n_kpts + 1,
                ),
                data,
            ]
        ),
        header=(
            "kpoint_index  abs_floquet_coefficient\n"
            f"eta={args.eta}, "
            f"initial_band={args.initial_band}, "
            f"projected_band={args.projected_band}"
        ),
    )

    print(
        f"[OK] Wrote coefficient data: "
        f"{coeff_file}"
    )

    if args.plot_kdist:
        plot_2D_kdist(
            data,
            lat,
            nspin=-1,
            plt_cbar=True,
            shift_BZ=False,
        )

        kdist_file = (
            f"floquet_kdist_eta_{args.eta:+d}_"
            f"b{args.initial_band}_to_b"
            f"{args.projected_band}.pdf"
        )

        plt.savefig(
            kdist_file,
            bbox_inches="tight",
        )

        plt.close()

        print(
            f"[OK] Wrote 2D k-distribution plot: "
            f"{kdist_file}"
        )


if __name__ == "__main__":
    main()
