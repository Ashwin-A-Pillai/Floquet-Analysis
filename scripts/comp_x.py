import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Input files
# ============================================================
file_NL = "o.YPP-eps_along_E"
file_RT = "o-CHI_IP.eps_q1_ip"

# ============================================================
# Read data
# ============================================================
# NL columns:
# E[eV], Im/eps_d1, Im/eps_d2, Im/eps_d3,
# Re/eps_d1, Re/eps_d2, Re/eps_d3
data_NL = np.loadtxt(file_NL, comments="#")

E_NL  = data_NL[:, 0]
Im_NL = data_NL[:, 1]
Re_NL = data_NL[:, 4]

# RT columns:
# E[eV], Im(eps), Re(eps)
data_RT = np.loadtxt(file_RT, comments="#")

E_RT  = data_RT[:, 0]
Im_RT = data_RT[:, 1]
Re_RT = data_RT[:, 2]


# ============================================================
# Plot 1: Imaginary part
# ============================================================
plt.figure(figsize=(7, 5))

plt.plot(
    E_NL, Im_NL,
    linewidth=2,
    label="Im/eps_d1 [NL]"
)

plt.plot(
    E_RT, Im_RT,
    linewidth=2,
    linestyle="--",
    label="Im(eps) [RT]"
)

plt.xlabel("E [eV]", fontsize=14)
plt.ylabel("Im($\\epsilon$)", fontsize=14)

plt.xlim(1.0, 10.0)

plt.legend(fontsize=12)
plt.tick_params(axis="both", labelsize=12)
plt.tight_layout()

plt.savefig("Im_eps_NL_RT.png", dpi=300, bbox_inches="tight")
plt.savefig("Im_eps_NL_RT.pdf", bbox_inches="tight")

plt.close()


# ============================================================
# Plot 2: Real part
# ============================================================
plt.figure(figsize=(7, 5))

plt.plot(
    E_NL, Re_NL,
    linewidth=2,
    label="Re/eps_d1 [NL]"
)

plt.plot(
    E_RT, Re_RT,
    linewidth=2,
    linestyle="--",
    label="Re(eps) [RT]"
)

plt.xlabel("E [eV]", fontsize=14)
plt.ylabel("Re($\\epsilon$)", fontsize=14)

plt.xlim(1.0, 10.0)

plt.legend(fontsize=12)
plt.tick_params(axis="both", labelsize=12)
plt.tight_layout()

plt.savefig("Re_eps_NL_RT.png", dpi=300, bbox_inches="tight")
plt.savefig("Re_eps_NL_RT.pdf", bbox_inches="tight")

plt.close()