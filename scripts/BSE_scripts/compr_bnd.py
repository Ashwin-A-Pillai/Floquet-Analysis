#!/usr/bin/env python3
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_GW_bands(filename):
    ks = []
    energies = []
    hs_k = []
    hs_labels = []

    band_min = band_max = None
    n_GW_bands = None

    with open(filename, 'r') as f:
        for ln in f:
            if '|k|' in ln:
                # header line: find b1...b8 ? digits [1,2,...,8]
                nums = re.findall(r'\d+', ln)
                band_min = int(nums[0])
                band_max = int(nums[-1])
                n_GW_bands  = band_max - band_min + 1
                continue

            if ln.lstrip().startswith('#') or not ln.strip():
                continue

            parts = ln.split()
            k = float(parts[0])
            ks.append(k)

            # next n_GW_bands columns are the energies
            rowE = [float(x) for x in parts[1:1+n_GW_bands]]
            energies.append(rowE)

            # check for trailing [LABEL]
            last = parts[-1]
            if last.startswith('[') and last.endswith(']'):
                hs_k.append(k)
                hs_labels.append(last.strip('[]'))

    # transpose energies into one list per band
    GW_bands = list(map(list, zip(*energies)))   # now GW_bands[i][j] = i-th band at j-th k

    return ks, GW_bands, hs_k, hs_labels


def plot_and_save(ks, GW_bands, KS_bands, hs_k, hs_labels, outname):
    fig, ax = plt.subplots(figsize=(8,6))

    # GW bands: label only the first
    for i, band in enumerate(GW_bands):
        lab = "GW" if i == 0 else "_nolegend_"
        ax.plot(ks, band, color="blue", label=lab, lw=1.5)
        
    # KS bands: label only the first
    for i, band in enumerate(KS_bands):
        lab = "Kohn-Sham" if i == 0 else "_nolegend_"
        ax.plot(ks, band, color="red", label=lab, lw=1.5)

    # vertical lines + ticks for high-symmetry
    if hs_k:
        for x in hs_k:
            ax.axvline(x=x, color='k', ls='--', lw=0.8)
        ax.set_xticks(hs_k)
        ax.set_xticklabels(hs_labels, fontsize=12)

    ax.set_xlabel('k-point', fontsize=14)
    ax.set_ylabel('Energy (eV)', fontsize=14)
    ax.grid(True, ls=':', lw=0.5)
    
    ax.legend()
    plt.tight_layout()
    fig.savefig(outname, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    fname1 = sys.argv[1] if len(sys.argv)>1 else 'GW_bands.out'
    fname2 = sys.argv[2] if len(sys.argv)>2 else 'KS_bands.out'
    outname = sys.argv[3] if len(sys.argv)>3 else 'bandstructure.png'

    ks, GW_bands, hs_k, hs_labels = read_GW_bands(fname1)
    ks, KS_bands, hs_k, hs_labels = read_GW_bands(fname2)
    plot_and_save(ks, GW_bands, KS_bands, hs_k, hs_labels, outname)