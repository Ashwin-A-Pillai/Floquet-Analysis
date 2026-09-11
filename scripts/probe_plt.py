#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Input file
# ---------------------------------------------------------------------------
input_file = Path("o.YPP-X_probe_order_2")

if not input_file.is_file():
    raise FileNotFoundError(f"Input file not found: {input_file}")

# Yambo output files are normally whitespace-separated.
# Lines beginning with "#" are ignored.
data = np.loadtxt(input_file, comments="#")

if data.ndim == 1:
    data = data.reshape(1, -1)

if data.shape[1] < 7:
    raise ValueError(
        f"Expected at least 7 columns, but found {data.shape[1]} "
        f"in {input_file}"
    )

# Gnuplot equivalent:
# u 1:(sqrt($6**2 + $7**2))
energy = data[:, 0]
chi2_abs = np.sqrt(data[:, 5] ** 2 + data[:, 6] ** 2)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 5.0))

ax.plot(
    energy,
    chi2_abs,
    linestyle="-",
    marker="o",
    markersize=5.5,
    linewidth=1.5,
)

ax.set_xlabel(r"Energy (eV)", fontsize=14)
ax.set_ylabel(
    r"$\left|\chi^{(2)}_{zxy}\right|\;(\mathrm{pm/V})$",
    fontsize=14,
)

ax.tick_params(axis="both", which="major", labelsize=12)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig("chi2_zxy.png", dpi=300, bbox_inches="tight")
fig.savefig("chi2_zxy.pdf", bbox_inches="tight")

plt.show()