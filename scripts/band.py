from collections import defaultdict
import re
import matplotlib.pyplot as plt 

#read file and get data
def read_GW(filename):

    file = open(filename,"r")
    data = defaultdict(list)

    band_max = 0
    band_min = 0

    for ln in file.readlines():

        chars = ln.split()

        if "|k|" in chars:
            
            nums = re.findall(r'\d+',ln)
            band_min = int(nums[0])
            band_max = int(nums[-1])
            continue
            
        if "#" not in chars:
            data["pkpt"].append(float(chars[0]))
            
            for i in range(1,band_max-band_min+2):
                data[f"band_{i+(band_max-band_min)-1}"].append(float(chars[i]))
        
        else:
            continue
    
    return data

def plot_bands(ax,data,color,linestyle="-"):
    for band in data:
        if band != "pkpt":
            ax.plot(data["pkpt"],data[band],linewidth=1.5,color=color,linestyle=linestyle)
            
            
# === main block ===
if __name__ == "__main__":
    # 1) create a figure & axes
    fig, ax = plt.subplots()
    
    fname = "bands[GW].out"
    pltname = "bands-GW.png"

    # 2) load and plot
    data = read_GW(fname)
    plot_bands(ax, data, color="blue", linestyle="-")

    # (optional) labels
    ax.set_xlabel("k-point")
    ax.set_ylabel("Energy (eV)")

    # 3) save to file instead of showing on screen
    fig.savefig(pltname, dpi=300, bbox_inches="tight")
    plt.close(fig)