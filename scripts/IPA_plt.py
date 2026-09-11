#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Input file
# ---------------------------------------------------------------------------
input_file = Path("CHI_IP/o-CHI_IP.eps_q1_ip")

if not input_file.is_file():
    raise FileNotFoundError(f"Input file not found: {input_file}")

# Columns:
# 1 -> Energy [eV]
# 2 -> Im(eps)
# 3 -> Re(eps)
data = np.loadtxt(input_file, comments="#")

energy = data[:, 0]
eps_imag = data[:, 1]
eps_real = data[:, 2]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 5.0))

ax.plot(
    energy,
    eps_imag,
    color="blue",
    linewidth=1.5,
    label=r"$\mathrm{Im}(\varepsilon)$",
)

ax.plot(
    energy,
    eps_real,
    color="red",
    linewidth=1.5,
    label=r"$\mathrm{Re}(\varepsilon)$",
)

ax.set_xlabel("Energy (eV)", fontsize=14)
ax.set_ylabel(r"$\varepsilon(\omega)$", fontsize=14)

ax.tick_params(axis="both", which="major", labelsize=12)

ax.legend(
    fontsize=12,
    frameon=False,
)

ax.grid(True, alpha=0.3)

fig.tight_layout()

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig.savefig("epsilon_IP.png", dpi=300, bbox_inches="tight")
fig.savefig("epsilon_IP.pdf", bbox_inches="tight")

plt.show()