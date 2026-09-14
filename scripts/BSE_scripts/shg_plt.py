#!/usr/bin/env python3

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Safe for SLURM / headless HPC nodes

import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------------
# Command-line argument
# -----------------------------------------------------------------------------
if len(sys.argv) != 2:
    print(f"Usage: python3 {sys.argv[0]} <YPP probe file>")
    sys.exit(1)

input_file = Path(sys.argv[1])

if not input_file.is_file():
    print(f"ERROR: File not found: {input_file}")
    sys.exit(1)


# -----------------------------------------------------------------------------
# Read Yambo/YPP output
# -----------------------------------------------------------------------------
# np.loadtxt automatically ignores lines beginning with '#'
data = np.loadtxt(input_file, comments="#")


if data.ndim == 1:
    data = data.reshape(1, -1)

if data.shape[1] < 5:
    print(
        f"ERROR: Expected at least 5 columns, "
        f"but found {data.shape[1]} in {input_file}"
    )
    sys.exit(1)


# -----------------------------------------------------------------------------
# Equivalent of:
#
# p 'o.YPP-X_probe_order_2' u 1:(sqrt($4**2+$5**2))
# -----------------------------------------------------------------------------
x = data[:, 0]

real_part = data[:, 3]   # Gnuplot column $4
imag_part = data[:, 4]   # Gnuplot column $5

shg = np.sqrt(real_part**2 + imag_part**2)


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(
    x,
    shg,
    marker="o",
    markersize=3,
    linewidth=1.5,
    label="SHG in hBN monolayer",
)

ax.set_xlabel("Energy")
ax.set_ylabel(r"$\sqrt{(\mathrm{Re})^2 + (\mathrm{Im})^2}$")

ax.legend()
ax.grid(True, alpha=0.3)

fig.tight_layout()


# -----------------------------------------------------------------------------
# Output file
# -----------------------------------------------------------------------------
output_file = input_file.parent / f"{input_file.name}_SHG.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)

print(f"SHG plot saved to: {output_file}")