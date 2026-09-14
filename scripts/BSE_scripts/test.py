import os, glob
from pathlib import Path

def head(path, n=40):
    try:
        with open(path, "r", errors="ignore") as f:
            for i in range(n):
                line = f.readline()
                if not line:
                    break
                print(line.rstrip())
    except Exception as e:
        print(f"  [could not read {path}: {e}]")

print("PWD:", os.getcwd())
print("\n=== Files here (top-level) ===")
for fn in sorted(os.listdir(".")):
    if any(fn.startswith(x) for x in ("ndb.", "o.", "r_", "s_", "l_")):
        print(" ", fn)

print("\n=== Looking for NL netCDF DBs (ndb.*) ===")
ndbs = sorted(glob.glob("ndb.*")) + sorted(glob.glob("**/ndb.*", recursive=True))
ndbs = [p for p in ndbs if os.path.isfile(p)]
print("Found", len(ndbs), "ndb.* files")
for p in ndbs[:50]:
    print(" ", p)
if len(ndbs) > 50:
    print(" ... (truncated)")

print("\n=== Looking for expected text outputs ===")
for p in ["o.external_field", "o.polarization", "r_nloptics"]:
    if os.path.exists(p):
        print("FOUND", p)
        if p.startswith("r_"):
            print("\n--- head:", p, "---")
            head(p, n=60)
    else:
        print("MISSING", p)

print("\n=== netCDF variable dump (if netCDF4 is available) ===")
try:
    import netCDF4
except Exception as e:
    print("netCDF4 import failed:", e)
    print("Install inside your conda env:  conda install -c conda-forge netcdf4")
    raise SystemExit(0)

def dump_one(ncpath, max_vars=200):
    print("\n---", ncpath, "---")
    ds = netCDF4.Dataset(ncpath, "r")
    try:
        print("Dimensions:", list(ds.dimensions.keys()))
        vars_ = list(ds.variables.keys())
        print("Num variables:", len(vars_))
        # show whether the key variable exists
        print("Has Field_Freq_range_1 ?", "Field_Freq_range_1" in ds.variables)
        print("Variables (first {}):".format(min(len(vars_), max_vars)))
        for v in vars_[:max_vars]:
            vv = ds.variables[v]
            dims = getattr(vv, "dimensions", ())
            shape = getattr(vv, "shape", ())
            print(f"  {v:35s} dims={dims} shape={shape}")
    finally:
        ds.close()

# Dump all top-level ndb.* first; if none, dump recursive matches
for p in sorted(glob.glob("ndb.*")):
    dump_one(p)
if not glob.glob("ndb.*"):
    for p in ndbs[:10]:
        dump_one(p)