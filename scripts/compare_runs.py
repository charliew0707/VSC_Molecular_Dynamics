"""
Compare runs pairwise: does a new run reproduce an old one?

For each pair it prints any in.json differences (ignoring absolute paths), how
far the x-dipole trajectories drift apart, and the two strongest IR peaks of
each (zero-padded FFT of the x dipole, the component `infrared -m ase` uses).

Usage:
    python scripts/compare_runs.py                    # the 2026-11-10 code check
    python scripts/compare_runs.py OLD_RUN NEW_RUN    # any two run folders
"""

import gzip
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
FS_TO_CM = 33356.41  # 1/(fs) in cm^-1

CODE_CHECK = [
    ("runs/co2/2026-10-02_lambda_sweep/nep_cav_w2400_lam0.1_chi",
     "runs/co2/2026-11-10_code_check/nep_cav_w2400_lam0.1_chi"),
    ("runs/co2/2026-10-02_lambda_sweep/nep_cav_w2400_lam0.1_nochi",
     "runs/co2/2026-11-10_code_check/nep_cav_w2400_lam0.1_nochi"),
    ("runs/h2o/2026-10-02_lambda_sweep/pyscf_cav_w1584_lam0.1_nochi",
     "runs/h2o/2026-11-10_code_check/pyscf_cav_w1584_lam0.1_nochi"),
    ("runs/h2o/2026-10-02_lambda_sweep/pyscf_cav_w1584_lam0.1_chi",
     "runs/h2o/2026-11-10_code_check/pyscf_cav_w1584_lam0.1_chi"),
]


def dipole(run):
    path = run / "dipole.dat"
    opener = open if path.exists() else gzip.open
    with opener(path if path.exists() else run / "dipole.dat.gz", "rt") as f:
        data = np.loadtxt(f)
    return data[:, 1], data[:, 2]  # time (a.u.), mu_x


def peaks(t, mu, lo=1000, hi=3500, n=2):
    dt = (t[1] - t[0]) * 0.0241888432658  # time column is in a.u.; convert to fs
    x = mu - mu.mean()
    spec = np.abs(np.fft.rfft(x, n=8 * len(x)))
    freq = np.fft.rfftfreq(8 * len(x), d=dt) * FS_TO_CM
    keep = (freq > lo) & (freq < hi)
    f, s = freq[keep], spec[keep]
    local = np.where((s[1:-1] > s[:-2]) & (s[1:-1] > s[2:]))[0] + 1
    top = local[np.argsort(s[local])[::-1][:n]]
    return sorted((round(f[i]), s[i] / s[top].max()) for i in top)


def inputs(run):
    with open(run / "in.json") as f:
        j = json.load(f)
    return {k: v for k, v in j.items() if not (isinstance(v, str) and v.startswith("/"))}


def compare(old, new):
    old, new = REPO / old, REPO / new
    print(f"\n{old.relative_to(REPO)}\n  vs {new.relative_to(REPO)}")
    a, b = inputs(old), inputs(new)
    diffs = {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)}
    print("  in.json differences:", diffs or "none")
    t1, m1 = dipole(old)
    t2, m2 = dipole(new)
    n = min(len(m1), len(m2))
    rel = np.abs(m1[:n] - m2[:n]).max() / np.abs(m1[:n]).max()
    print(f"  steps: {len(m1) - 1} vs {len(m2) - 1};  max |d mu_x| / max |mu_x| = {rel:.2e}")
    hi = 2200 if "h2o" in str(old) else 3500
    for label, (t, m) in (("old", (t1, m1)), ("new", (t2, m2))):
        print(f"  {label} peaks (cm^-1, relative height):",
              ", ".join(f"{f} ({h:.2f})" for f, h in peaks(t, m, hi=hi)))


if __name__ == "__main__":
    pairs = [tuple(sys.argv[1:3])] if len(sys.argv) == 3 else CODE_CHECK
    for old, new in pairs:
        compare(old, new)
