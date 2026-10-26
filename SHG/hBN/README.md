# Nonlinear optics and Floquet analysis of monolayer hBN

**Quantum ESPRESSO → Lumen/Yambo → YamboPy** · Independent-particle and screened-exchange calculations · Second-harmonic generation (SHG) · Floquet–Kohn–Sham states · Excitons

This directory contains the inputs, analysis programs, selected numerical datasets and figures from two complementary calculations of **monolayer hexagonal boron nitride (hBN)**:

1. **IPA** — independent-particle approximation, with no explicit time-dependent electron–hole interaction.
2. **TD-SHF / TD-HSEX** — time-dependent screened Hartree–Fock (equivalently, in the present implementation, time-dependent Hartree plus screened exchange), including the screened electron–hole interaction and a quasiparticle scissor correction.

The purpose is to compare the **real-time (RT) nonlinear polarization** with a **Floquet reconstruction** of the same driven dynamics, examine the resulting quasienergies and Fourier coefficients, and correlate excitonic levels and momentum-space weights with the observed SHG resonances. The published numerical arrays can be loaded directly by third-party Python scripts: the supplied PNGs are not the only means of visualization.

> **Scientific scope.** An RT/Floquet comparison tests the internal consistency of the reconstruction; it is not, by itself, a convergence study in k points, KS bands, driving intensity or Floquet truncation. The reported nonlinear coefficients are *finite-field, field-contracted effective responses* unless a separately identifiable tensor element is selected. Proximity between a BSE exciton energy and an SHG peak does **not** prove that exciton's individual contribution to SHG.

## Contents

- [Calculations and physical parameters](#calculations-and-physical-parameters)
- [Directory layout](#directory-layout)
- [The three Slurm jobs](#the-three-slurm-jobs)
- [Input files and Python tools](#input-files-and-python-tools)
- [Numerical data dictionary](#numerical-data-dictionary)
- [Figure guide](#figure-guide)
- [Make your own plots](#make-your-own-plots)
- [Interpretation and reproducibility](#interpretation-and-reproducibility)
- [References](#references)

## Calculations and physical parameters

| Setting | [IPA](hBN/IPA/) | [TD-SHF / TD-HSEX](hBN/TDSHF/) |
|:--|:--|:--|
| Material | Monolayer hBN (2 atoms per cell) | Monolayer hBN (same structural model) |
| Electronic-structure code | Quantum ESPRESSO | Quantum ESPRESSO |
| SCF/NSCF k-point mesh | **30 × 30 × 1** | **18 × 18 × 1** |
| Plane-wave kinetic-energy cutoff | **80 Ry** | **80 Ry** |
| NSCF KS bands (`nbnd`) | **200** | **200** |
| Real-time evolution | Lumen `yambo_nl` | Lumen `yambo_nl` |
| Correlation option | `NLCorrelation="IPA"` | `NLCorrelation="SEX"` (Hartree + screened exchange) |
| Quasiparticle correction in archived RT/BSE inputs | `GfnQP_E=0` | **3.64670676 eV** scissor (`GfnQP_E` for RT, `KfnQP_E` for BSE) |
| Screening/BSE/collision preparation | Not used for IPA evolution | `GW0.in`, `BSE.in`, `COL.in` |
| Propagated active occupied KS states | Bands 3–4 | Bands 3–4 |
| KS basis for propagation | `NLBands=3..6` | `NLBands=3..6` |
| Driving field | Linear, along Cartesian **y**, `SOFTSIN` | Linear, along Cartesian **y**, `SOFTSIN` |
| Field intensity as recorded by Lumen | `Field1_Int=1000.0 kWLm2` | `Field1_Int=1000.0 kWLm2` |
| Time step | `NLstep=0.01 fs` | `NLstep=0.01 fs` |
| Input dephasing parameter | `NLDamping=0.20 eV` | `NLDamping=0.20 eV` |
| Real-time driving energies | **191** nodes, approximately **0.5–10.0 eV** | **191** nodes, approximately **0.5–10.0 eV** |
| Selected Floquet driving energies | **64** adaptive RT-grid nodes | **64** adaptive RT-grid nodes |
| Stored Fourier/Floquet order | **`FLOrder=4`** | **`FLOrder=4`** |
| Main observable | SHG, the `y` output for a `y`-polarized field | Same, plus exciton-resonance analysis |

The RT frequencies and the 64 Floquet frequencies are **different sampling grids of driving photon energy**; they are not electronic band energies. At each selected driving energy, Lumen propagates the Bloch states for that field frequency. A Floquet order of four yields **nine Fourier sidebands** with indices `η = −4, …, 4`. This value was not increased to conceal numerical errors.

The two calculations have **different k meshes** and **different Hamiltonian approximations**. They are not a controlled single-parameter comparison of screened exchange alone. For the TD-SHF/BSE data, the quoted scissor is a **conduction-band correction**, not the exciton energy; the first bright BSE excitonic doublet is near **6.00 eV**.

## Directory layout

Paths below are relative to this `SHG/README.md` file.

```text
SHG/
├── README.md                         # This guide
├── FLQ_AL.py                        # Common database, lattice, Berry, tensor and audit utilities
├── FLQuasi.py                       # Floquet eigenproblem + Fourier reconstruction
├── FLQ_HG.py                        # RT and reconstructed Floquet harmonic susceptibilities
├── FLQPlot.py                       # Spectrum, Floquet/KS, time-series and BZ plots
├── FLQ_Exc.py                       # BSE exciton catalogue, resonance matches and k-space maps
├── patch_*.sh                       # Input patchers for Lumen/Yambo and YPP
├── hBN/
│   ├── IPA/
│   │   ├── Inputs/                  # QE, RT, Floquet and YPP example inputs
│   │   ├── SLURM/                   # Three published HPC job scripts
│   │   ├── Data/
│   │   │   ├── PyFLQ/               # Complete spectra, comparisons, logs
│   │   │   ├── FLQ_023/PyFLQ/       # One representative frequency: Floquet state data
│   │   │   └── o-DE2.YPP-X_probe_order_2
│   │   └── *.png                    # Prepared spectrum/state illustrations
│   └── TDSHF/
│       ├── Inputs/                  # QE, GW0, BSE, COL, RT, Floquet and YPP inputs
│       ├── Data/
│       │   ├── PyFLQ/               # Complete spectra, comparisons, logs
│       │   ├── FLQ_026/PyFLQ/       # One representative frequency: Floquet state data
│       │   ├── ExcFLQ/              # BSE exciton weights, tables and summaries
│       │   └── o-DE2.YPP-X_probe_order_2
│       └── *.png                    # SHG, Floquet and exciton figures
├── Paper/                            # Draft manuscript, figures and references
└── AlAs/                             # Separate material; outside this hBN guide
```

**Published-data limitation:** the repository contains the *full RT/Floquet SHG spectra*, but only **one example per-frequency Floquet state archive for each method** (`FLQ_023` for IPA and `FLQ_026` for TD-SHF). It does **not** presently include all 64 large per-frequency raw Lumen databases, the full `SAVE/ns.db1`, the BSE eigenvector database, or a complete self-contained HPC runtime environment. The high-symmetry-path CSV referenced by the original scripts (`../../../IBZ/hBN.csv` relative to the run directory) is also **not present in the public repository snapshot**; the path coordinates can instead be obtained from the saved `bands_path.npz`.

## The three Slurm jobs

The **three committed `.job` files are under [`hBN/IPA/SLURM/`](hBN/IPA/SLURM/)**. They implement a staged workflow, with a separate submission for each stage. TD-SHF uses an analogous three-stage RT → frequency-array → postprocessing workflow, but **the TD-SHF-specific Slurm scripts are not committed under `hBN/TDSHF/`**. Its actual `GW0.in`, `BSE.in`, `COL.in`, `YNL.in` and sample Floquet inputs are provided there.

| Job script | Role | Key products |
|:--|:--|:--|
| [`NL_final.job`](hBN/IPA/SLURM/NL_final.job) | **Setup + full RT calculation.** Converts the QE NSCF database with `p2y`, sets up `SAVE`, calculates linear IPA optics, runs YPP `FixSymm`, performs the RT nonlinear frequency sweep (`NL`), extracts the RT harmonic spectrum, selects an adaptive 64-point Floquet frequency grid, and prepares the Floquet template. | `SAVE`, `FixSymm`, `FS_IPA`, `YNL.in`, `NL`, `PyFLQ/RT_order_2.dat`, `FLQ_frequencies.dat`, `FLQ_frequencies.json`, `.flq.env` |
| [`FrqArr.job`](hBN/IPA/SLURM/FrqArr.job) | **Parallel frequency-resolved Lumen Floquet propagation.** Reads the prepared 64-frequency list, generates `Y_FLQ_###.in` with `patch_FP.sh`, and starts independent 8-rank `yambo_nl` runs with `srun`. Uses a Slurm array and concurrent slots within each task. | `FLQ_000/`, …, `FLQ_063/`, Lumen `ndb.RT_V_bands*`, `ndb.Nonlinear`, `.complete` marker |
| [`YPY.job`](hBN/IPA/SLURM/YPY.job) | **Postprocessing and plotting.** Checks the completion signatures, runs YPP harmonic analysis, extracts and audits Floquet quasienergies/Fourier modes, reconstructs Berry-phase SHG, compares against RT SHG, plots both spectra and one selected Floquet state/BZ example. | `PyFLQ/HG_order_2.npz`, `Floquet_order_2.dat`, plots, and per-frequency `quasienergies.npz`, `fourier.npz`, `HG_order_2.npz` |

**Important execution notes:**

- These are HPC scripts using local CINaM partitions, environment modules, absolute executable paths, and a Conda environment. Adapt paths and modules before running elsewhere.
- `NL_final.job` explicitly deletes/rebuilds selected setup directories (including `SAVE`, `FixSymm`, `CHI_IP`) during fresh preparation. Do **not** rerun it on irreplaceable original datasets without a backup.
- `FrqArr.job` assigns **four frequencies per array task** (`--ntasks=32`, 8 ranks per frequency) and limits the Slurm array to five simultaneous tasks (`--array=0-15%5`). A `.complete` marker causes a frequency directory to be skipped; an existing directory without it is deleted and rerun. Never change the frequency list while reusing the old indexed `FLQ_###` folders.
- `YPY.job` executes its scientific Python analysis loops in series in this committed version. Merely increasing CPU allocation will not make those loops parallel. The script typically produces only **one example state/bandstructure plot**, chosen by `--target-eV` (default 4.1 eV).
- `.flq_branch`, `.flq.env`, `.RT_input.sha256` and `.complete` are **run-provenance and orchestration files**, not scientific observables.

## Input files and Python tools

### Quantum ESPRESSO and Lumen/Yambo inputs

| File(s) | Description |
|:--|:--|
| `Inputs/hBN.scf.K*.E80.B200.in` | Quantum ESPRESSO self-consistent electronic-structure setup. |
| `Inputs/hBN.nscf.K*.E80.B200.in` | NSCF calculation supplying wavefunctions, 200 bands, 80 Ry and the specified 30×30×1 or 18×18×1 mesh. |
| `IPA/Inputs/IPA.in` | Linear-optics calculation in the independent-particle approximation (not the Floquet pump calculation). |
| `TDSHF/Inputs/GW0.in` | GW₀ stage and screening-related setup; existence of GW₀ does **not** mean that its `ndb.QP` corrections were used in the released scissor-mode RT results. |
| `TDSHF/Inputs/BSE.in` | Bethe–Salpeter equation, screened electron–hole kernel, exciton energies and oscillator strengths; the `KfnQP_E` scissor is ~3.6467 eV. |
| `TDSHF/Inputs/COL.in` | Screened-exchange/collision-kernel preparation needed by the interacting RT simulation. |
| `Inputs/FSypp.in` | YPP symmetry-reduction (`FixSymm`) input: preserves operations compatible with the driven electric field. Do not use this reduced symmetry group to classify the **equilibrium** SHG tensor. |
| `Inputs/YNL.in` | Full multi-frequency nonlinear RT calculation; `NLEnRange`, `NLEnSteps`, `NLBands`, `NLCorrelation`, damping and field settings. |
| `Inputs/Y_FLQ_039.in` | Example **single-frequency** Floquet propagation input (`FLOrder=4`). The `039` index is a sample input and is **not** the archived `FLQ_023`/`FLQ_026` state sample. |
| `Inputs/anl.in` | YPP harmonic fitting/postprocessing input. |

The shared root-level [`patch_*.sh`](.) programs modify autogenerated Yambo/Lumen inputs without performing a simulation: `patch_IPA.sh` (linear IPA), `patch_GW.sh` (GW), `patch_BSE.sh` (BSE), `patch_FSy.sh` (symmetry fixing), `patch_COL.sh` (collision/screened exchange), `patch_FN.sh` (RT `NL` sweep), `patch_FP.sh` (single-frequency Floquet), and `patch_FAL.sh` (YPP harmonic analysis). They expect specific environment variables; consult the Slurm scripts and patchers before reproducing the workflow.

### Analysis scripts

| Script | What it does |
|:--|:--|
| [`FLQ_AL.py`](FLQ_AL.py) | Central utilities: Lumen/NetCDF and `.npz` I/O, k-mesh and symmetry handling, Berry-phase polarization, tensor-component classification, RT/Floquet matching, numerical tests, extraction audit, and adaptive Floquet frequency selection. |
| [`FLQuasi.py`](FLQuasi.py) | Builds the one-period monodromy matrix for the **propagated occupied subspace**, obtains quasienergies and Floquet eigenvectors, checks one-/four-period closure, and saves Fourier/Sambe amplitudes. |
| [`FLQ_HG.py`](FLQ_HG.py) | Computes RT harmonic susceptibilities and reconstructs Floquet-state Berry-phase polarization; performs a temporal sampling check and builds matched RT-versus-Floquet SHG datasets. |
| [`FLQPlot.py`](FLQPlot.py) | Generates susceptibility curves, complex real/imaginary panels, RT/Floquet differences, Floquet-versus-KS band plots, individual wavefunction projections, and BZ Fourier-weight maps. |
| [`FLQ_Exc.py`](FLQ_Exc.py) | Uses BSE energies, oscillator strengths and eigenvectors to identify one-photon-weak states, plot normalized k-space exciton weights, and mark **candidate** one-/two-photon resonances in SHG. |

## Numerical data dictionary

### A. Full RT and Floquet SHG spectra — `Data/PyFLQ/`

**Start with these files for a custom SHG plot**:

- [`IPA/Data/PyFLQ/HG_order_2.npz`](hBN/IPA/Data/PyFLQ/HG_order_2.npz)
- [`TDSHF/Data/PyFLQ/HG_order_2.npz`](hBN/TDSHF/Data/PyFLQ/HG_order_2.npz)

The main arrays have the same schema in both methods:

| Key | Shape / type | Meaning |
|:--|:--|:--|
| `rt_energy_eV` | `(191,)` | Pump photon energy of each full real-time calculation, in eV. |
| `rt_chi` | `(191, 3)`, complex | Effective harmonic response at each RT energy; output Cartesian components **x, y, z** in that order. |
| `fl_energy_eV` | `(64,)` | Selected driving energies for the Floquet calculation, in eV. |
| `fl_chi` | `(64, 3)`, complex | Same complex effective harmonic response, reconstructed from the saved Floquet states. |
| `rt_matched_chi` | `(64, 3)`, complex | RT values selected at the corresponding Floquet energies for point-by-point comparison. |
| `rt_matching_indices` | `(64,)`, integer | Zero-based indices into `rt_energy_eV` that generated `rt_matched_chi`. |
| `rt_selected_chi`, `fl_selected_chi` | `(191,)`, `(64,)`, complex | A selected observable (here the `y`-output field contraction), possibly a symmetry-identifiable tensor element in other runs. |
| `rt_raw_harmonics` | `(fit_order+1, 191, 3)`, complex | RT temporal Fourier-fit coefficients prior to response normalization. |
| `harmonic_order`, `fit_order` | scalar integers | Requested optical harmonic (here 2) versus Fourier-fit cutoff (here 4). **Not `FLOrder`.** |
| `selected_component`, `selected_label`, `selected_kind` | strings | Selection and identification, e.g. `y`, `chi_y_eff`, `field_contracted`. |
| `units` | string | For SHG in the released files, `cm/statV`. |
| `normalization`, `thickness_factor` | string, scalar | Supercell-volume convention and any optional slab-thickness rescaling. Here spectra report unrescaled supercell response. |
| `tensor_catalog_json` | JSON string | Equilibrium symmetry/tensor classification and which tensor elements a single field orientation can identify. |
| `metadata_json` | JSON string | Input parameters, per-frequency diagnostics, RT/Floquet comparison metrics and provenance information. |

**Do not compare every RT sample with the Floquet arrays by array index.** The Floquet grid is an **adaptive subset**, with only 64 points; use `rt_matched_chi` or `rt_matching_indices` to compare matching energies. `fl_chi[:, 1]` and `rt_chi[:, 1]` are the emitted **y** components for the incident y-polarized field. This is generally written as `chi_y_eff`; it should not be relabelled as an independently computed full Cartesian tensor.

The same response is published in human-readable form: `RT_order_2.dat`, `Floquet_order_2.dat` (seven numerical columns: `E_eV`, `Im_x`, `Re_x`, `Im_y`, `Re_y`, `Im_z`, `Re_z`), and their `*_selected.dat` versions (three columns: `Energy_eV`, `Im_response`, `Re_response`). `Data/o-DE2.YPP-X_probe_order_2` is the original YPP-generated harmonic reference output. `*.log` files report checks, settings, extraction status and errors; `plots/*_metrics.json` stores quantitative RT/Floquet comparison measures.

### B. Per-frequency Floquet-state cache — `Data/FLQ_###/PyFLQ/`

Published representative directories:

- [IPA `FLQ_023`](hBN/IPA/Data/FLQ_023/PyFLQ/)
- [TD-SHF `FLQ_026`](hBN/TDSHF/Data/FLQ_026/PyFLQ/)

`###` is a zero-based *position in the selected frequency list*, **not an energy in eV**. The sample archived objects use different field frequencies and are not expected to correspond across the two methods by index.

**`quasienergies.npz`** — Floquet eigenproblem and equilibrium electronic reference:

| Key | Dimensions | Interpretation |
|:--|:--|:--|
| `epsilon_principal_au`, `xi_principal_eV`, `multipliers` | `(Nk_IBZ, 2)` | Principal quasienergies (atomic units/eV) and one-period Floquet multipliers. Quasienergies are periodic modulo `ħω`. |
| `L`, `L4`, `Gamma`, `Gram` | `(Nk_IBZ, 2, 2)` | One-period and four-period projected evolution, Floquet eigenvectors, and initial-state overlap metric. |
| `mode_for_band`, `mode_assignment_weight` | `(Nk_IBZ, 2)` | Mapping from the two propagated occupied KS bands to the two represented Floquet modes and assignment confidence. Mode indices here are **zero-based**. |
| `kpoints_red`, `kpoints_full_red` | `(Nk_IBZ, 3)`, `(Nk_BZ, 3)` | Reduced-coordinate k-points in the **field-reduced** irreducible zone and full periodic mesh. |
| `bz_to_ibz` | `(Nk_BZ,)`, integer | Zero-based full-BZ → IBZ map. |
| `ks_energies_eV` | `(Nk_IBZ, 200)` | KS band energies, in eV, on the driven IBZ mesh. |
| `ks_vbm_eV` | scalar | KS valence-band maximum used as a plot reference. |
| `grid_shape`, `varying_axes`, `lattice_bohr`, `reciprocal_no2pi_bohr` | arrays | k-grid topology, lattice geometry and reciprocal-space conventions. |
| `omega_au`, `period_au`, `t0_au`, `nf` | scalars | Driving angular frequency, period, reference time and saved Floquet order (`nf=4`). |
| `active_bands`, `basis_bands` | `(2,)`, `(4,)` | **One-based** KS indices: active = `[3,4]`; basis = `[3,4,5,6]`. |
| `closure_T`, `closure_4T`, `semigroup_4T`, `gram_unitarity`, `gauge_*`, ... | `(Nk_IBZ,)` | Numerical quality checks; see [`FLQuasi.py`](FLQuasi.py) for their definitions. |
| `schema`, `generation`, `metadata_json` | strings | Cache schema, generation UUID and validation/provenance. |

The archived per-frequency examples have `Nk_IBZ=465, Nk_BZ=900` (IPA) and `Nk_IBZ=171, Nk_BZ=324` (TD-SHF). These are **field-reduced k meshes**; they differ from the **equilibrium** BSE IBZ used in the exciton maps.

**`fourier.npz`** — periodic Floquet components and original checkpoints:

| Key | Shape | Interpretation |
|:--|:--|:--|
| `Lambda` | `(Nk_IBZ, 2, 9, 4)`, complex | Floquet–KS amplitudes `Λ[k, mode, eta, basis_band]`. The nine `eta` entries correspond to `−4..4`. |
| `harmonics` | `(9,)` | Actual harmonic indices `[-4,-3,-2,-1,0,1,2,3,4]`; do not assume the array index equals `eta`. |
| `epsilon_compact_au`, `gauge_shift_m` | `(Nk_IBZ,2)` | Quasienergies in a compact Fourier gauge and their integer branch shifts. |
| `edge_fraction`, `spectral_moment` | `(Nk_IBZ,2)` | Diagnostics of weight near the truncation boundary and sideband distribution. |
| `C_raw` | `(Ntime, Nk_IBZ, 4, 2)`, complex | Saved propagated KS coefficients at the Lumen checkpoint times. |
| `raw_times_au`, `raw_indices`, `phase_sample_indices`, `stored_reconstruction_errors` | 1D arrays | Original time/checkpoint bookkeeping and Fourier reconstruction errors. |
| `omega_au`, `period_au`, `t0_au`, `generation` | shared scalars/identifier | Must agree with `quasienergies.npz`. |

The state convention in the script is schematically `ψ_a(t0+τ) = exp(−i ε_a τ) Σ_η Λ[a,η] exp(+i η ω τ)` in atomic units, followed by a transformation with the inverse `Gamma` matrix to reconstruct the propagated occupied orbitals. **Two reconstructed occupied modes do not constitute a complete Floquet conduction-plus-valence bandstructure.** KS-connected integer quasienergy branches and global reference alignment used for display are visualization conventions, not QP corrections.

**`HG_order_2.npz` (inside one `FLQ_###`)** contains the **single-frequency** reconstructed response: `energy_eV` (scalar), `chi` (three complex Cartesian components), `raw_harmonics` (`(fit_order+1,3)`), `times_au`, `polarization` (`(Ntime,3)`), `projector` (`(3,3)`), and `metadata_json`. It is not the same as the aggregated `Data/PyFLQ/HG_order_2.npz`.

**`HG_order_2_Berry_debug.npz`** records `times_au`, `raw_unprojected`, `raw_projected`, `stored` and `projector`. The two polarization arrays are reconstructed from saved Floquet coefficients and compared with the original Lumen polarization; use them to study Berry reconstruction and symmetry leakage, not as an exciton-wavefunction database.

### C. Arrays associated with the state and Brillouin-zone plots

| Filename | Main keys and shapes | Use |
|:--|:--|:--|
| `FLQ_###/PyFLQ/plots/bands_path.npz` | `distance_inv_ang (Npath,)`, `points_red (Npath,3)`, `ks_energies_eV (Npath,8)`, `energies_eV (Npath,2)`, `ks_energies_relative_to_reference_eV`, `floquet_energies_relative_to_reference_eV`, tick labels/positions and alignment metadata | Replot the **two occupied reconstructed Floquet modes** (plus optional ±replicas) with all eight KS reference bands along Γ–K–M–Γ. Published paths have `Npath=91`. |
| `FLQ_###/PyFLQ/plots/FKS_BZ.npz` | `kpoints_full_cart_inv_ang (Nk_BZ,3)`, `values_full (Nk_BZ,)`, `bz_vertices`, `plot_frame`, `quantity`, `selected_modes_1based`, `selected_bands`, `bz_to_ibz` | Map the selected Fourier/Sambe weight in the first BZ. For the default plotted case, it is a sum of `|Λ|²` at `η=+2` over represented modes and KS basis bands, **not an exciton probability**. |
| `FLQ_###/PyFLQ/plots/plot_settings.json` | Human-readable plotting arguments | Recreate the figure settings: KS overlay, Floquet gauge, global phase-reference alignment, band range, harmonic `η`, and high-symmetry path choice. |

The full complex `Lambda` array from `fourier.npz`, rather than `FKS_BZ.npz`, is required for a new band-, mode-, or harmonic-resolved Fourier analysis. Mode/band-resolved full-BZ maps may need symmetry sewing information near degeneracies; the committed plotting code defaults to complete sums for that reason.

### D. Excitons and dark-state candidates — `TDSHF/Data/ExcFLQ/`

The archived BSE analysis covers **5,184 excitons**. The first bright doublet is at about **6.000 eV**; states 5, 6 and 9 are **weak in the calculated one-photon BSE channel**, not automatically two-photon-allowed states.

| Files | Description |
|:--|:--|
| [`exciton_catalogue.csv`](hBN/TDSHF/Data/ExcFLQ/exciton_catalogue.csv) | One row per exciton: one-based state index, BSE energy, normalized one-photon strength, candidate pump energies `E_exc/2` and `E_exc`, and weak-state flag. |
| [`multiplet_resonances.csv`](hBN/TDSHF/Data/ExcFLQ/multiplet_resonances.csv) | Selected nearly degenerate exciton groups, mean energies/splittings, summed strengths and candidate resonance energies. |
| [`SHG_local_peaks.csv`](hBN/TDSHF/Data/ExcFLQ/SHG_local_peaks.csv) | Numerically detected RT-SHG magnitude peaks and complex susceptibility values. |
| [`SHG_exciton_resonance_matches.csv`](hBN/TDSHF/Data/ExcFLQ/SHG_exciton_resonance_matches.csv) | Nearest sampled SHG values/peaks to BSE `E_exc/2` and `E_exc`; **energetic coincidence, not an exciton decomposition**. |
| [`dark_candidates_SHG.csv`](hBN/TDSHF/Data/ExcFLQ/dark_candidates_SHG.csv) | Candidate one-photon-weak excitons and their nearest SHG features. |
| [`dark_vs_bright_resolvability.csv`](hBN/TDSHF/Data/ExcFLQ/dark_vs_bright_resolvability.csv) | Bright–weak resonance separation relative to a **user-supplied broadening scale**, not a measured linewidth. |
| [`analysis_summary.json`](hBN/TDSHF/Data/ExcFLQ/analysis_summary.json) | Input paths, selected exciton groups, peak summaries, numerical status, and explicit interpretation warnings. |
| `weights_XXXX[_YYYY].npz` | Normalized per-exciton or averaged-nearly-degenerate-subspace momentum-space weights. See below. |
| `weights_*_IBZ.csv` | Irreducible-star positions, star multiplicity, integrated weight and per-star-point weight. |
| `weights_*_HS_path.csv` | Γ–K–M–Γ sampled path: path distance and **nearest computed mesh point** and its weight. It is not a continuously computed exciton dispersion. |
| `weights_*_transitions.csv` | Band-pair-resolved mean probabilities, labelled by one-based valence and conduction KS band. |

**`weights_0005.npz` or `weights_0001_0002.npz`** provides:

| Key | Shape for this 18×18×1 hBN BSE | Meaning |
|:--|:--|:--|
| `states` | `(1,)` or `(2,)` | One-based BSE exciton IDs represented by the file. |
| `IBZ_kred` | `(37,3)` | Equilibrium BSE **irreducible** reduced-coordinate k points. |
| `IBZ_star_size` | `(37,)` | Number of full-BZ points in each equilibrium symmetry star. |
| `IBZ_star_weight` | `(37,)` | **Integrated probability per star**; sums to 1 for each saved group-averaged map. |
| `full_kred` | `(324,3)` | Full 18×18×1 BSE mesh, in reduced coordinates. |
| `full_kweight` | `(324,)` | **Probability per full-BZ k point**; sums to 1. |

For a nearly degenerate group of excitons `G`, the maps use `W_G(k) = (1/|G|) Σ_{λ∈G} Σ_{vc} |A_{λ,vck}|²`. Averaging the full group prevents arbitrary rotations within a nearly degenerate subspace from changing the reported distribution. `IBZ_star_weight` and `full_kweight` are **not interchangeable**: the first integrates over a symmetry star; the second is normalized **per full-zone point**. Neither array stores the complex phase of the BSE eigenvectors. The published data therefore support exciton **probability maps**, but not independent reconstruction of coherent two-photon or exciton–exciton optical matrix elements without the original BSE/dipole databases.

## Figure guide

### IPA and TD-SHF SHG spectra

Each method publishes these figures directly in `hBN/<method>/`:

| Image | What is plotted / how to interpret it |
|:--|:--|
| `HG_order_2_chi_y_eff_magnitude.png` | `|χ^(2)_eff,y|` versus **driving** photon energy: continuous RT spectrum and adaptively sampled Floquet points. Units `cm/statV`. |
| `HG_order_2_chi_y_eff_complex.png` | **Top:** `Re χ^(2)_eff,y`; **bottom:** `Im χ^(2)_eff,y`. Sign changes and dispersive/resonant structures are visible here even when magnitude alone is ambiguous. |
| `HG_order_2_chi_y_eff_difference.png` | RT–Floquet **magnitude difference** normalized to the peak RT magnitude, in percent. This is not the same as a complex relative error. |
| `FLQ_adaptive_frequency_grid.png` | Distribution of the **64** Floquet driving energies chosen from the **191** RT nodes. The additional panel is an estimated interpolation mismatch used in choosing nodes, not a simulated SHG error bar. |

[IPA magnitude](hBN/IPA/HG_order_2_chi_y_eff_magnitude.png) · [IPA complex](hBN/IPA/HG_order_2_chi_y_eff_complex.png) · [TD-SHF magnitude](hBN/TDSHF/HG_order_2_chi_y_eff_magnitude.png) · [TD-SHF complex](hBN/TDSHF/HG_order_2_chi_y_eff_complex.png)

![TD-SHF: RT and Floquet SHG magnitude](hBN/TDSHF/HG_order_2_chi_y_eff_magnitude.png)

### Representative Floquet/KS and Fourier-state plots

The following images exist in **both** `hBN/IPA/` and `hBN/TDSHF/` and correspond to the representative frequency saved in each `Data/FLQ_###/PyFLQ/` directory:

| Image | Interpretation |
|:--|:--|
| `fig-Floquet_bands_IBZ.png` | Floquet quasienergies of the **two represented occupied modes**, with integer `mħω` replicas, alongside **KS bands 1–8** on Γ–K–M–Γ. The energy zero is the KS VBM; the visualization employs KS-connected gauge and a documented **global phase-reference alignment**. A Floquet replica is not a distinct independently propagated conduction state. |
| `fig-FKS_projection.png` | Floquet–Kohn–Sham coefficient magnitudes `|Λ|` for a selected k point, initially associated occupied band, KS basis band and Fourier sideband `η`. |
| `fig-FKS-proj_BZ.png` | Full-first-BZ distribution of a selected **Floquet sideband weight** (typically summed `|Λ|²` at `η=2`); it is *not* the BSE exciton wavefunction. |
| `fig-real_time_projection_Re.png` / `_Im.png` | Real and imaginary parts of a demodulated periodic Floquet-mode coefficient: reconstructed curve versus stored Lumen checkpoints. |
| `fig-propagated_state_Re.png` / `_Im.png` | Real and imaginary parts of an **original propagated occupied state** reconstructed from the full mode combination, versus stored checkpoints. |

[IPA band figure](hBN/IPA/fig-Floquet_bands_IBZ.png) · [TD-SHF band figure](hBN/TDSHF/fig-Floquet_bands_IBZ.png) · [IPA BZ sideband map](hBN/IPA/fig-FKS-proj_BZ.png) · [TD-SHF BZ sideband map](hBN/TDSHF/fig-FKS-proj_BZ.png)

![TD-SHF Floquet modes and KS reference](hBN/TDSHF/fig-Floquet_bands_IBZ.png)

### TD-SHF exciton and dark-state plots

All files below are published under [`hBN/TDSHF/`](hBN/TDSHF/). `0001_0002` denotes **BSE excitons 1 and 2** grouped together (one-based IDs); `0005` denotes exciton 5 only.

| Figure | What the image communicates |
|:--|:--|
| [`01_BSE_exciton_strengths.png`](hBN/TDSHF/01_BSE_exciton_strengths.png) | BSE eigenenergies and normalized **one-photon** oscillator strengths; near-zero values identify candidates dark for the calculated polarization. |
| [`02_SHG_BSE_resonances.png`](hBN/TDSHF/02_SHG_BSE_resonances.png) | RT SHG data with one-photon (`E_exc`) and two-photon (`E_exc/2`) *candidate resonance markers*. Marker coincidence does **not** assign the nonlinear amplitude. |
| [`02_SHG_resonance_guide.png`](hBN/TDSHF/02_SHG_resonance_guide.png) | Separate schematic/numerical reference of BSE energies and potential driving-energy resonances. |
| [`03_dark_exciton_SHG_zoom.png`](hBN/TDSHF/03_dark_exciton_SHG_zoom.png) | Zoom of the SHG spectral region near bright and one-photon-weak candidate excitons, including states 5 and 6. Nearness to the broad peak is **not a resolved dark-exciton signature**. |
| `04_IBZ_states_*.png` | **37-point equilibrium IBZ** plots of **symmetry-star-integrated** exciton weights, one group per file. |
| `05_BZ_states_*.png` | Spatial distribution of normalized **full-BZ** electron–hole transition probability over the 18×18×1 mesh, one exciton group per file. |
| `06_HS_path_states_*.png` | Excitonic weight sampled along **Γ–K–M–Γ** using the nearest available BSE k point; an intensity profile, **not** an exciton dispersion `E_exc(k)`. |
| [`07_dark_vs_bright_IBZ.png`](hBN/TDSHF/07_dark_vs_bright_IBZ.png) | IBZ comparison of a bright multiplet (1,2) and a one-photon-weak state (5), and their weight difference. This highlights changes in momentum-space structure but does not reveal eigenvector phases. |

The available `04`, `05` and `06` families cover **`0001_0002`, `0003_0004`, `0005`, `0006`, `0007_0008`, `0009`, `0010_0011`**. Example links: [bright BZ map](hBN/TDSHF/05_BZ_states_0001_0002.png), [weak-state 5 BZ map](hBN/TDSHF/05_BZ_states_0005.png), [state 5 IBZ map](hBN/TDSHF/04_IBZ_states_0005.png), [state 5 HS path](hBN/TDSHF/06_HS_path_states_0005.png).

![Bright and weak exciton weights in the BSE IBZ](hBN/TDSHF/07_dark_vs_bright_IBZ.png)

## Make your own plots

A local Python environment needs **NumPy** and **Matplotlib** for the following examples. None of them requires installing Lumen or YamboPy because they read **published `.npz` and CSV arrays** only. Run them from the `SHG/` directory. The arrays are saved without Python object pickles, so use `allow_pickle=False`.

### Example 1 — SHG magnitude and complex response from the complete dataset

```python
import numpy as np
import matplotlib.pyplot as plt

file = "hBN/TDSHF/Data/PyFLQ/HG_order_2.npz"  # or hBN/IPA/...
with np.load(file, allow_pickle=False) as data:
    E_rt = data["rt_energy_eV"]      # (191,) in eV
    E_fl = data["fl_energy_eV"]      # (64,) in eV
    chi_rt = data["rt_chi"][:, 1]    # y output, complex
    chi_fl = data["fl_chi"][:, 1]
    units = str(data["units"])

fig, ax = plt.subplots()
ax.plot(E_rt, np.abs(chi_rt), label="RT")
ax.scatter(E_fl, np.abs(chi_fl), s=14, label="Floquet")
ax.set(xlabel="Driving photon energy (eV)",
       ylabel=rf"$|\chi^{{(2)}}_{{\mathrm{{eff}},y}}|$ ({units})")
ax.legend()
fig.tight_layout()
plt.show()
```

For a complex response, use `chi_rt.real` and `chi_rt.imag`. To get the matched comparison rather than interpolate manually, use `data["rt_matched_chi"][:, 1]` with `data["fl_chi"][:, 1]`.

### Example 2 — KS bands and reconstructed Floquet modes along Γ–K–M–Γ

```python
import numpy as np
import matplotlib.pyplot as plt

path = "hBN/TDSHF/Data/FLQ_026/PyFLQ/plots/bands_path.npz"
with np.load(path, allow_pickle=False) as d:
    x = d["distance_inv_ang"]
    ks = d["ks_energies_relative_to_reference_eV"]
    floq = d["floquet_energies_relative_to_reference_eV"]
    tick_x = d["tick_positions"]
    tick_names = d["tick_labels"]

fig, ax = plt.subplots()
for b in range(ks.shape[1]):
    ax.plot(x, ks[:, b], color="tab:red", lw=1, label="KS" if b == 0 else None)
for a in range(floq.shape[1]):
    ax.plot(x, floq[:, a], "--", color="tab:blue", lw=1.2,
            label="Occupied Floquet modes" if a == 0 else None)
ax.set_xticks(tick_x, tick_names)
ax.set(xlabel="High-symmetry path", ylabel="Energy relative to KS VBM (eV)")
ax.legend()
fig.tight_layout()
plt.show()
```

This draws **the central Floquet representatives only**. To recreate the replica curves in the publication image, add `m * ℏω` for integer `m` (here from `−2` to `+2`) using `omega_au × 27.211386245988` from the matching `quasienergies.npz`. Do not interpret the copies as extra computed eigenstates.

### Example 3 — Replot a BSE exciton weight map in reduced k coordinates

```python
import numpy as np
import matplotlib.pyplot as plt

file = "hBN/TDSHF/Data/ExcFLQ/weights_0005.npz"
with np.load(file, allow_pickle=False) as d:
    k = d["full_kred"]           # (324,3), reduced reciprocal coordinates
    weight = d["full_kweight"]   # (324,), sum = 1
    states = d["states"]         # [5], one-based BSE exciton ID

fig, ax = plt.subplots()
points = ax.scatter(k[:, 0], k[:, 1], c=weight, s=35)
fig.colorbar(points, ax=ax, label="Probability per full-BZ k point")
ax.set(xlabel="Reduced k1", ylabel="Reduced k2",
       title=f"Exciton {states.tolist()}: transition weight")
ax.set_aspect("equal")
fig.tight_layout()
plt.show()
```

This last example shows the **periodic full mesh in reduced coordinates**; it does not apply Wigner–Seitz first-BZ clipping or Cartesian k-vector conversion. See [`FLQ_Exc.py`](FLQ_Exc.py) and `05_BZ_states_*.png` for the first-BZ presentation. The `weights` files are **probability distributions**; they cannot reconstruct phase winding or exciton-to-exciton dipoles without the missing complex BSE eigenvectors.

### Discover every `.npz` key / shape yourself

```python
from pathlib import Path
import numpy as np

for file in Path("hBN").rglob("*.npz"):
    print("\n", file)
    with np.load(file, allow_pickle=False) as db:
        for key in db.files:
            a = db[key]
            print(f"  {key:42} shape={a.shape!s:16} dtype={a.dtype}")
```

## Interpretation and reproducibility

**RT versus Floquet SHG.** The output is the same finite-field effective second-order response evaluated in two ways: directly from the real-time polarization and by reconstructing the periodic Floquet states and their Berry polarization. The data include both **complex susceptibility** and **magnitude**; phases matter near resonances. The published `HG_order_2_chi_y_eff_metrics.json` gives complex L² relative discrepancies of approximately **0.50% (IPA)** and **0.70% (TD-SHF)** at the 64 matched frequencies. These are cross-reconstruction metrics, not basis/k-mesh convergence bars.

**Excitonic SHG.** The TD-SHF data have a lowest bright BSE exciton doublet at approximately **6.0 eV** and one-photon-weak candidates near **6.94 eV** (state 5) and **7.00 eV** (state 6). The full RT susceptibility shows prominent magnitude features near pump energies **3.00, 3.45 and 5.95 eV**. Conditions `2ℏω≈E_exc` and `ℏω≈E_exc` are useful *resonance guides*. A strong SHG peak close to a one-photon-dark exciton does **not** demonstrate that the exciton is two-photon active or uniquely contributes to the peak; that requires symmetry/selection-rule analysis and coherent interexciton optical matrix elements. The IBZ/BZ maps plot the probability `Σ_vc|A^λ_vck|²`, not an exciton-resolved `χ^(2)`.

**Symmetry.** `SAVE` before the field-induced symmetry reduction is the correct source for **equilibrium** tensor selection rules; `FixSymm/SAVE` is appropriate to the **field-compatible real-time propagation** and to expanding driven-state coefficients. In the archived hBN point group, the symmetry-allowed second-order tensor subspace has one independent parameter. Do not infer this from the reduced group alone.

**Bands and energies.** KS energies and Floquet quasienergies are different objects. Quasienergies are defined **modulo the pump photon energy** and have gauge freedom; the KS-connected/global-offset choices used for the figures are plotting conventions. The archived `Lambda` contains only **two dynamically represented occupied modes**, with frozen deeper occupied bands treated separately in the Berry determinant. This is not an independently computed full valence-and-conduction Floquet eigenspectrum.

**What is not contained in this public snapshot:** original electronic `SAVE/ns.db1`, the full collection of raw `ndb.RT_V_bands*` at all frequencies, original `ndb.BS_diago_Q1` eigenvectors/`ndb.dipoles`, and the local `IBZ/hBN.csv`. Consequently a collaborator **can replot the released data** without Lumen, but **cannot rerun all stages or calculate new coherent interexciton amplitudes solely from this repository**. The `Paper/` directory contains a draft manuscript and background references, not necessarily an accepted or final publication.

## References

- [Quantum ESPRESSO](https://www.quantum-espresso.org/) — SCF/NSCF electronic structure.
- [Lumen](https://www.lumen-code.org/) and [Yambo](https://www.yambo-code.eu/) — GW, BSE, screening, real-time nonlinear optics and Floquet-state sampling.
- [YamboPy Floquet-analysis tutorial (`fl-analysis` branch)](https://github.com/yambo-code/yambopy/blob/fl-analysis/yambopy/nl/vb_postprocessing_tutorial/Floquet_analysis_tutorial.md) — example real-time → Floquet workflow.
- [YamboPy exciton/BZ postprocessing tutorial](https://www.lumen-code.org/wiki/index.php?title=Exciton-phonon_coupling_and_luminescence_-_Yambopy_postprocessing) and [How to analyse excitons](https://www.lumen-code.org/wiki/index.php?title=How_to_analyse_excitons).
- [Theory and work-in-progress manuscript](Paper/main.tex) — direct monodromy/Floquet–KS formulation used for the projected, actively propagated occupied subspace.

If you use the published arrays, please cite the computational packages and any associated eventual project publication when one becomes available; no finished-paper citation is implied by this repository snapshot.
