#!/usr/bin/env python3
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_bands(filename):
    ks = []
    energies = []
    hs_k = []
    hs_labels = []

    band_min = band_max = None
    n_bands = None

    with open(filename, 'r') as f:
        for ln in f:
            if '|k|' in ln:
                # header line: find b1...b8 ? digits [1,2,...,8]
                nums = re.findall(r'\d+', ln)
                band_min = int(nums[0])
                band_max = int(nums[-1])
                n_bands  = band_max - band_min + 1
                continue

            if ln.lstrip().startswith('#') or not ln.strip():
                continue

            parts = ln.split()
            k = float(parts[0])
            ks.append(k)

            # next n_bands columns are the energies
            rowE = [float(x) for x in parts[1:1+n_bands]]
            energies.append(rowE)

            # check for trailing [LABEL]
            last = parts[-1]
            if last.startswith('[') and last.endswith(']'):
                hs_k.append(k)
                hs_labels.append(last.strip('[]'))

    # transpose energies into one list per band
    bands = list(map(list, zip(*energies)))   # now bands[i][j] = i-th band at j-th k

    return ks, bands, hs_k, hs_labels


def plot_and_save(ks, bands, hs_k, hs_labels, outname):
    fig, ax = plt.subplots(figsize=(8,6))

    # plot each band
    for band in bands:
        ax.plot(ks, band, color="blue", lw=1.5)

    # vertical lines + ticks for high-symmetry
    if hs_k:
        for x in hs_k:
            ax.axvline(x=x, color='k', ls='--', lw=0.8)
        ax.set_xticks(hs_k)
        ax.set_xticklabels(hs_labels, fontsize=12)

    ax.set_xlabel('k-point', fontsize=14)
    ax.set_ylabel('Energy (eV)', fontsize=14)
    ax.grid(True, ls=':', lw=0.5)

    plt.tight_layout()
    fig.savefig(outname, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv)>1 else 'bands.out'
    outname = sys.argv[2] if len(sys.argv)>2 else 'bandstructure.png'

    ks, bands, hs_k, hs_labels = read_bands(fname)
    plot_and_save(ks, bands, hs_k, hs_labels, outname)