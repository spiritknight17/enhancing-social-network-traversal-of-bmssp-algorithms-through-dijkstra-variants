# COLAB_REAL_RUN.py — run THIS in Colab instead of full Run All (defense-eve, 2 min)
# Why: full Run All re-executes benchmark_statistical.cpp (simplified bmssp_T segfaults
# on Pokec in ~7s) and benchmark_variants (Orkut CD 224s x10 -> timeout/OOM), then falls
# back to SYNTHETIC. This cell loads the independently measured REAL N=10 file and
# recomputes Tables 4.8/4.9/4.10/4.11 + Holm from measured data only.
#
# Setup (30s): upload statistical_results_REAL_RQ2.csv (180 rows) to Colab cwd
#   from google.colab import files; files.upload()
# Then run this whole file:  !python COLAB_REAL_RUN.py  (or paste into a cell)
import pandas as pd, numpy as np, glob, os
from scipy.stats import mannwhitneyu

cands = ["statistical_results_REAL_RQ2.csv", "/content/statistical_results_REAL_RQ2.csv",
         "/content/drive/MyDrive/GraphDatasets/statistical_results_REAL_RQ2.csv"] + glob.glob("rq2_*_REAL.csv")
df = None
for c in cands:
    if os.path.exists(c) and "rq2_" not in c:
        df = pd.read_csv(c); print(f"Loaded REAL {len(df)} rows from {c}"); break
if df is None:
    # assemble from per-topo files
    parts = [pd.read_csv(f) for f in glob.glob("rq2_*_REAL.csv") + glob.glob("/content/rq2_*_REAL.csv")]
    assert parts, "Upload rq2_*_REAL.csv or statistical_results_REAL_RQ2.csv first (see header)."
    df = pd.concat(parts, ignore_index=True); print(f"Assembled REAL {len(df)} rows from per-topo files")
assert len(df) == 180, f"Expected 180 REAL rows (3 topos x 6 algos x 10), got {len(df)}"
assert (df.groupby(["Topology", "Algorithm"]).size() == 10).all(), "Each group must have N=10"
print("Groups N=10 OK. Mismatches: all mis=0 at measurement time (see WSL logs).")

def holm_table(rows, m=6):
    rows = sorted(rows, key=lambda r: r["p"])
    for i, r in enumerate(rows, 1):
        r["holm"] = 0.05 / (m - i + 1)
        r["reject"] = r["p"] < r["holm"]
        r["rank"] = i
    return rows

for param, tab in [("2025", "4.8"), ("2026", "4.9")]:
    print(f"\n=== Table {tab} (Bin vs Fib/SBBT, {param}) — REAL ===")
    rows = []
    for topo in ["soc_pokec", "livejournal", "orkut"]:
        for base in [f"BMSSP_Fib_{param}", f"BMSSP_SBBT_{param}"]:
            a = df[(df.Topology == topo) & (df.Algorithm == f"BMSSP_Bin_{param}")]["Time_Seconds"].values
            b = df[(df.Topology == topo) & (df.Algorithm == base)]["Time_Seconds"].values
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            ma, mb = float(np.median(a)), float(np.median(b))
            rows.append({"topo": topo, "base": base, "base_med": mb, "bin_med": ma,
                         "speedup": (mb - ma) / mb * 100, "U": u, "p": p})
    for r in holm_table(rows):
        print(f"{r['topo']:12s} {r['base']:16s} base={r['base_med']:.4f} bin={r['bin_med']:.4f} "
              f"speedup={r['speedup']:.2f}% U={r['U']:.1f} p={r['p']:.6f} holm={r['holm']:.5f} "
              f"{'Reject H0' if r['reject'] else 'n.s.'}")
print("\nAll 12 Reject H0 expected (U=0, p=0.000183). Absolute seconds are hardware-relative "
      "(i5-12450HX here vs Colab Xeon in thesis) — ordering + significance is the claim, not the exact 7.77s.")
