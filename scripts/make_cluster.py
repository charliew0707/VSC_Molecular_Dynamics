"""
Build a small cluster of N identical, aligned molecules for multi-molecule runs.

Copies a single-molecule geometry N times, shifting copy k by k*R along an axis, and
writes an extended-XYZ file into geometries/ with a name that says what it is:

    geometries/<mol>_N<N>_<arrangement>-<axis>_R<R>A.xyz
    e.g. geometries/h2o_N2_headtail-x_R8.0A.xyz

Arrangements (cavity polarization and each molecule's dipole are along x):
    headtail  copies shifted along x, the dipole axis (molecules in a line, end to end)
    sidebyside copies shifted along y, perpendicular to the dipoles

R is the distance between neighbouring molecules' centres (same shift for every atom,
so every copy has exactly the starting geometry and orientation: only the bright
collective mode is excited). The box is only a label for cboamd (non-periodic runs),
sized so the cluster fits with at least 10 A to spare on each side.

Usage:
    python scripts/make_cluster.py h2o 2 8.0                     # head-to-tail along x
    python scripts/make_cluster.py h2o 2 4 6 8 10                # several separations
    python scripts/make_cluster.py h2o 3 8.0 --arrangement sidebyside
"""

import argparse
import json
from pathlib import Path

import numpy as np
from ase import io as aio

REPO = Path(__file__).resolve().parent.parent
AXES = {"headtail": (0, "x"), "sidebyside": (1, "y")}


def build(molecule, n, r, arrangement):
    cfg = json.loads((REPO / "configs" / "molecules" / f"{molecule}.json").read_text())
    source = REPO / cfg["xyz_file"]
    single = aio.read(source)
    axis, axis_name = AXES[arrangement]

    shift = np.zeros(3)
    shift[axis] = r
    cluster = single.copy()
    for k in range(1, n):
        copy = single.copy()
        copy.translate(k * shift)
        cluster += copy
    cluster.translate(-cluster.get_positions().mean(axis=0))  # centre the cluster at 0

    extent = cluster.get_positions().max(axis=0) - cluster.get_positions().min(axis=0)
    box = max(20.0, float(np.ceil(extent.max() + 20.0)))
    cluster.set_cell([box, box, box])
    cluster.set_pbc(False)
    cluster.info = {"molecules": n, "molecule": molecule, "arrangement": arrangement,
                    "axis": axis_name, "R_angstrom": r, "source": source.name,
                    "made_by": "scripts/make_cluster.py"}

    name = f"{molecule}_N{n}_{arrangement}-{axis_name}_R{r:.1f}A.xyz"
    out = REPO / "geometries" / name
    if out.exists():
        print(f"exists, left unchanged: geometries/{name}")
        return out
    aio.write(out, cluster, format="extxyz")
    print(f"wrote geometries/{name}  ({len(cluster)} atoms, box {box:.0f} A)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("molecule", help="a configs/molecules/<name>.json (its xyz_file is copied)")
    ap.add_argument("n", type=int, help="number of molecules")
    ap.add_argument("r", type=float, nargs="+", help="separation(s) R in angstrom")
    ap.add_argument("--arrangement", choices=AXES, default="headtail")
    args = ap.parse_args()
    for r in args.r:
        build(args.molecule, args.n, r, args.arrangement)


if __name__ == "__main__":
    main()
