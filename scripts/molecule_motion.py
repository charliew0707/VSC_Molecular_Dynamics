"""
How much each molecule rotates and drifts during a run.

For every molecule in the run (the molecule count comes from the "molecules=N" entry
that scripts/make_cluster.py writes into the xyz header; plain single-molecule files
count as 1) it reports, over the whole trajectory:

    rot   largest angle (deg) the molecule's axis turns away from where it started.
          Axis = first atom -> centre of the others (for H2O: O -> midpoint of the
          H atoms, i.e. the C2 axis, along which the bend dipole points). Molecules
          in a tilted field swing back and forth, so this is the largest swing,
          not the final angle.
    tilt  largest angle (deg) between that axis and the cavity polarization
          (lambda_vector; x when there is no cavity).
    drift largest displacement (A) of the molecule's centre of mass.
and for clusters the smallest and largest centre-to-centre distance (A).

Reads position.dat (written every step by cboamd), else md.traj. Both stay on the
machine that ran the job (not in git), so this only works there.

Usage:
    python scripts/molecule_motion.py runs/h2o_N2/<campaign>          # every run in it
    python scripts/molecule_motion.py runs/h2o_N2/<campaign>/<run> ...
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase import io as aio

EVERY = 5   # analyse every 5th step (2.5 fs at dt = 0.5 fs): plenty for rotations


def trajectory(run, natoms):
    """Positions, shape (frames, natoms, 3), from position.dat or md.traj."""
    if (run / "position.dat").exists():
        data = np.loadtxt(run / "position.dat")[::EVERY]
        return data[:, 2:2 + 3 * natoms].reshape(len(data), natoms, 3)
    if (run / "md.traj").exists():
        return np.array([a.positions for a in aio.read(run / "md.traj", f"::{EVERY}")])
    return None


def motion(run):
    cfg = json.loads((run / "in.json").read_text())
    start = aio.read(cfg["xyz_file"])
    nmol = int(start.info.get("molecules", 1))
    natoms = len(start)
    per = natoms // nmol
    pos = trajectory(run, natoms)
    if pos is None:
        return None
    masses = start.get_masses()
    pol = np.asarray(cfg.get("lambda_vector", [[1, 0, 0]])[0], dtype=float)
    pol /= np.linalg.norm(pol)

    out = dict(molecules=nmol, frames=len(pos), mol=[])
    centres = []
    for m in range(nmol):
        sl = slice(m * per, (m + 1) * per)
        p, w = pos[:, sl], masses[sl]
        axis = p[:, 1:].mean(axis=1) - p[:, 0]
        axis /= np.linalg.norm(axis, axis=1, keepdims=True)
        rot = np.degrees(np.arccos(np.clip(axis @ axis[0], -1, 1))).max()
        tilt = np.degrees(np.arccos(np.clip(np.abs(axis @ pol), 0, 1))).max()
        com = (p * w[None, :, None]).sum(axis=1) / w.sum()
        drift = np.linalg.norm(com - com[0], axis=1).max()
        centres.append(com)
        out["mol"].append(dict(rot=rot, tilt=tilt, drift=drift))
    if nmol > 1:
        c = np.array(centres)                       # (molecules, frames, 3)
        d = [np.linalg.norm(c[i] - c[j], axis=1) for i in range(nmol) for j in range(i + 1, nmol)]
        out["sep_min"], out["sep_max"] = float(np.min(d)), float(np.max(d))
    out["max_rot"] = max(x["rot"] for x in out["mol"])
    out["max_drift"] = max(x["drift"] for x in out["mol"])
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    runs = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        runs += [p] if (p / "in.json").exists() else sorted(d for d in p.iterdir() if (d / "in.json").exists())
    print(f"{'run':52s} {'mol':>3s} {'rot':>6s} {'tilt':>6s} {'drift':>6s} {'sep min-max':>13s}")
    for run in runs:
        m = motion(run)
        if m is None:
            print(f"{run.name:52s}  (no position.dat or md.traj here)")
            continue
        for i, x in enumerate(m["mol"]):
            sep = f"{m['sep_min']:5.2f}-{m['sep_max']:5.2f}" if (i == 0 and "sep_min" in m) else ""
            print(f"{run.name if i == 0 else '':52s} {i + 1:3d} {x['rot']:6.1f} {x['tilt']:6.1f} "
                  f"{x['drift']:6.3f} {sep:>13s}")
    print("\nrot/tilt in degrees, drift and separation in angstrom")


if __name__ == "__main__":
    main()
