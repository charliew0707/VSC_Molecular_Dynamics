"""
Energy-conservation check for cboamd runs (e-mode, one cavity mode).

The conserved quantity is
    H = Ekin(nuclei) + E_bare + 1/2 p^2 + 1/2 A ea^2,    A = 1 + lam^2 (pol.chi.pol)
where md.log's Etot = Ekin + E_bare (it does NOT include the cavity), p and ea are
the photon momentum and (screened) field from photon.dat, and A = 1 when the
polarizability term is off (polar=false). With A = 1 and chi on, H looks like it
fluctuates ~8%; that's a bookkeeping artefact, not the dynamics.

Reports, as fractions of the run's largest kinetic + photon energy:
    fluct  max |H - <H>|   (good: ~1%)
    drift  linear-fit change of H over the whole run   (good: < 0.5%)

Usage:
    python scripts/check_energy.py runs/h2o/2026-10-02_lambda_sweep [more run or campaign dirs]
Needs the full outputs (md.log, photon.dat, polarizability.dat), which stay on
the machine that ran them (not in git).
"""

import json
import sys
from pathlib import Path

import numpy as np

HA = 27.211386        # eV per hartree
FLUCT_OK, DRIFT_OK = 0.03, 0.005


def check(run):
    cfg = json.loads((run / "in.json").read_text())
    log = np.loadtxt(run / "md.log", skiprows=1)
    etot, ekin = log[:, 1], log[:, 3]
    n = len(etot)
    if cfg.get("photons"):
        ph = np.loadtxt(run / "photon.dat")
        n = min(n, len(ph))
        p, ea = ph[:n, 3], ph[:n, 4]
        a = np.ones(n)
        if cfg.get("polar", True):
            pol = np.loadtxt(run / "polarizability.dat")[:n, 2:11].reshape(-1, 3, 3)
            e = np.asarray(cfg.get("lambda_vector", [[1, 0, 0]])[0], dtype=float)
            e /= np.linalg.norm(e)
            a = 1 + cfg["lambda_photon"][0] ** 2 * np.einsum("i,nij,j->n", e, pol, e)
        eph = (0.5 * p ** 2 + 0.5 * a * ea ** 2) * HA
    else:
        eph = np.zeros(n)
    h = etot[:n] + eph
    scale = (ekin[:n] + eph - eph.min()).max()
    fluct = np.abs(h - h.mean()).max() / scale
    drift = np.polyfit(np.arange(n), h - h[0], 1)[0] * (n - 1) / scale
    return n - 1, fluct, drift


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    runs = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        runs += [p] if (p / "in.json").exists() else sorted(d for d in p.iterdir() if (d / "md.log").exists())
    print(f"{'run':52s} {'steps':>6s} {'fluct':>7s} {'drift':>8s}")
    bad = 0
    for run in runs:
        steps, fluct, drift = check(run)
        ok = fluct < FLUCT_OK and abs(drift) < DRIFT_OK
        bad += not ok
        print(f"{run.parent.name + '/' + run.name:52s} {steps:6d} {fluct:7.4f} {drift:+8.4f}  {'ok' if ok else 'CHECK'}")
    print(f"\n{len(runs) - bad}/{len(runs)} runs conserve energy (fluct < {FLUCT_OK:.0%}, |drift| < {DRIFT_OK:.1%})")


if __name__ == "__main__":
    main()
