"""
Plot a two-water separation scan (made by scripts/campaigns/h2o_N2_separation_scan.sh).

Question: from which separation R do two waters couple only through the cavity? Then
  - the bare bend (no cavity) is the same as for one water, and
  - the Rabi splitting is sqrt(2) times the single-molecule one at the same lambda
    (collective coupling of the bright mode).

Reads every run in the campaign folder: cavity runs (name ..._R<R>) and bare runs
(name ..._R<R>_bare); R comes from the geometry's xyz header. The single-molecule
reference is the 2026-10-02 H2O lambda sweep at the same lambda (same level of theory,
same 800-step analysis). Peaks use the same method as plot_h2o_sweep.py (8x zero-padded
FFT of dipole x, second peak must be >= 50% of the first). At 800 steps (83 cm^-1
resolution) peak positions are good to about +-10 cm^-1 (tested on synthetic data with
planted peaks), so compare ratios to sqrt(2) = 1.41 with that in mind.

Molecule motion (largest rotation and centre-of-mass drift, closest approach) comes
from scripts/molecule_motion.py, which needs the full outputs on this machine; it is
saved in peaks.csv so the plot can be remade from a fresh clone.

Usage:
    python scripts/plot_h2o_N2_separation.py runs/h2o_N2/2026-11-18_separation_scan
    python scripts/plot_h2o_N2_separation.py <campaign> --out <dir>   # write elsewhere (tests)

Writes to results/h2o_N2/<campaign>/:
    separation.png   (a) Rabi splitting vs R with 1 x and sqrt(2) x the single-molecule
                     value, (b) bare bend vs R with the single-molecule bend
    ir_spectra.png   stacked bend-region spectra per R, bare and in the cavity
    peaks.csv        every run: R, lambda, bend / LP / UP / Rabi, rotation, drift
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase import io as aio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import molecule_motion as mm                      # noqa: E402
import plot_h2o_sweep as sw                       # noqa: E402  (style, FFT, peak finder)

REPO = sw.REPO
REFERENCE = REPO / "runs" / "h2o" / "2026-10-02_lambda_sweep"
COLOR, MARKER = sw.COLOR, sw.MARKER               # PySCF, chi neglected: same series colour
INK, INK2 = sw.INK, sw.INK2
CM_PER_AU = 219474.63
FIELDS = ["run", "R_A", "lambda", "steps", "resolution_cm-1", "bend_cm-1", "LP_cm-1", "UP_cm-1",
          "Rabi_cm-1", "resolved", "max_rotation_deg", "max_drift_A", "min_separation_A"]


def dipole_file(d):
    return d / "dipole.dat" if (d / "dipole.dat").exists() else d / "dipole.dat.gz"


def analyse(d):
    cfg = json.loads((d / "in.json").read_text())
    dip = np.loadtxt(dipole_file(d))
    if len(dip) < cfg["steps"] + 1:
        print(f"skipping {d.name}: incomplete ({len(dip)}/{cfg['steps'] + 1} steps)")
        return None
    geom = REPO / "geometries" / Path(cfg["xyz_file"]).name   # in.json paths are absolute
    info = aio.read(geom).info
    freq, amp = sw.spectrum(dip[:, 2], cfg["timestep"])
    cavity = bool(cfg.get("photons"))
    lp, up, resolved = sw.polaritons(freq, amp)
    return dict(run=d.name, path=d, R=float(info.get("R_angstrom", np.nan)),
                molecules=int(info.get("molecules", 1)), cavity=cavity,
                lam=cfg["lambda_photon"][0] if cavity else 0.0, steps=cfg["steps"],
                res=1 / (cfg["steps"] * cfg["timestep"] * 1e-15) / sw.C_CM_S,
                omega_c=cfg["omega_photon"][0] * CM_PER_AU if cavity else None,
                freq=freq, amp=amp, bend=sw.local_maxima(freq, amp, sw.BEND_BAND)[0][0],
                lp=lp, up=up, resolved=resolved and cavity)


def reference(lam):
    """Single water at the same lambda (and bare), from the 2026-10-02 sweep."""
    runs = {r["lam"]: r for r in sw.load_campaign(REFERENCE, "nochi")}
    if lam not in runs:
        sys.exit(f"no single-molecule reference at lambda = {lam} in {REFERENCE}")
    return runs[0.0]["lp"], runs[lam]["up"] - runs[lam]["lp"]


def add_motion(runs, cached):
    for r in runs:
        m = mm.motion(r["path"])
        if m is not None:
            r.update(rot=m["max_rot"], drift=m["max_drift"], sep_min=m.get("sep_min", np.nan))
        else:
            c = cached.get(r["run"], {})
            r.update(rot=float(c.get("max_rotation_deg") or np.nan),
                     drift=float(c.get("max_drift_A") or np.nan),
                     sep_min=float(c.get("min_separation_A") or np.nan))


def plot_separation(cav, bare, bend1, rabi1, out):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    lam = cav[0]["lam"] if cav else None
    if cav:
        R = [r["R"] for r in cav]
        x0, x1 = min(R) - 0.5, max(R) + 0.3
        a.hlines(rabi1, x0, x1, color=INK2, ls=":", lw=1.3)
        a.hlines(np.sqrt(2) * rabi1, x0, x1, color=INK2, ls="--", lw=1.3)
        a.text(x1 + 0.1, rabi1, f"1 molecule\n({rabi1:.0f})", va="center", fontsize=8.5, color=INK2)
        a.text(x1 + 0.1, np.sqrt(2) * rabi1, f"√2 × 1 molecule\n({np.sqrt(2) * rabi1:.0f})",
               va="center", fontsize=8.5, color=INK2)
        ok = [r for r in cav if r["resolved"]]
        a.plot([r["R"] for r in ok], [r["up"] - r["lp"] for r in ok], ls="-", lw=1.2, color=COLOR,
               marker=MARKER, ms=9, label=f"2 H₂O, λ = {lam:g} (PySCF, χ neglected)")
        for r in cav:
            if not r["resolved"]:
                a.plot(r["R"], 0, marker="x", color=INK2, ms=8)
        a.axhspan(0, cav[0]["res"], color="#edece8", zorder=0)
        a.text(min(R), cav[0]["res"] * 0.45, f"below resolution ({cav[0]['res']:.0f} cm⁻¹)",
               fontsize=8, color=INK2)
        a.set_xlim(min(R) - 0.5, max(R) + 3.2)
        a.set_ylim(0, max(np.sqrt(2) * rabi1, max([r["up"] - r["lp"] for r in ok] or [0])) * 1.25)
    a.set_xlabel("Separation R (Å, centre to centre)")
    a.set_ylabel("Rabi splitting Ω$_R$ = UP − LP (cm⁻¹)")
    a.set_title("(a) Splitting vs separation", loc="left", fontsize=10)
    a.legend(loc="lower right", fontsize=8.5)

    if bare:
        R = [r["R"] for r in bare]
        b.hlines(bend1, min(R) - 0.5, max(R) + 0.3, color=INK2, ls=":", lw=1.3)
        b.text(max(R) + 0.4, bend1, f"1 molecule\n({bend1:.0f})", va="center", fontsize=8.5, color=INK2)
        b.plot(R, [r["bend"] for r in bare], ls="-", lw=1.2, color=INK, marker="o", ms=7,
               label="2 H₂O, no cavity")
        b.axhspan(bend1 - bare[0]["res"] / 2, bend1 + bare[0]["res"] / 2, color="#edece8", zorder=0)
        b.set_xlim(min(R) - 0.5, max(R) + 2.6)
        b.legend(loc="lower right", fontsize=8.5)
    b.set_xlabel("Separation R (Å, centre to centre)")
    b.set_ylabel("Bend peak (cm⁻¹)")
    b.set_title("(b) Bare bend vs separation (band = ± half the resolution)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "separation.png", dpi=150)
    plt.close(fig)


def plot_spectra(cav, bare, out):
    Rs = sorted({r["R"] for r in cav + bare})
    fig, ax = plt.subplots(figsize=(7, 1.2 + 1.1 * len(Rs)))
    for i, R in enumerate(Rs):
        for r, color, label in [(next((x for x in bare if x["R"] == R), None), INK2, "no cavity"),
                                (next((x for x in cav if x["R"] == R), None), COLOR, "cavity")]:
            if r is None:
                continue
            m = (r["freq"] > sw.BEND_BAND[0]) & (r["freq"] < sw.BEND_BAND[1])
            ax.plot(r["freq"][m], i + 0.85 * r["amp"][m] / r["amp"][m].max(), color=color, lw=1.6,
                    label=label if i == 0 else None)
        ax.text(sw.BEND_BAND[0] + 15, i + 0.55, f"R = {R:g} Å", fontsize=9)
    omega = next((r["omega_c"] for r in cav), None)
    if omega:
        ax.axvline(omega, color=INK2, ls="--", lw=1)
        ax.text(omega + 8, len(Rs) - 0.15, "ω$_c$", fontsize=9, color=INK2)
    ax.set_yticks([])
    ax.set_xlabel("Frequency (cm⁻¹)")
    ax.set_ylabel("IR intensity (normalized, offset by R)")
    ax.set_title("Two H₂O, bend region: no cavity vs cavity", loc="left", fontsize=10)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out / "ir_spectra.png", dpi=150)
    plt.close(fig)


def main():
    args = sys.argv[1:]
    if len(args) not in (1, 3) or (len(args) == 3 and args[1] != "--out"):
        sys.exit(__doc__)
    campaign = Path(args[0]).resolve()
    out = Path(args[2]).resolve() if len(args) == 3 else REPO / "results" / campaign.parent.name / campaign.name
    out.mkdir(parents=True, exist_ok=True)
    cached = {}
    if (out / "peaks.csv").exists():
        with open(out / "peaks.csv") as f:
            cached = {row["run"]: row for row in csv.DictReader(f)}

    runs = [r for d in sorted(campaign.iterdir()) if (d / "in.json").exists() and (r := analyse(d))]
    add_motion(runs, cached)
    cav = sorted((r for r in runs if r["cavity"]), key=lambda r: r["R"])
    bare = sorted((r for r in runs if not r["cavity"]), key=lambda r: r["R"])
    lams = {r["lam"] for r in cav}
    if len(lams) > 1:
        sys.exit(f"campaign mixes lambda values {sorted(lams)}; plot one lambda at a time")
    bend1, rabi1 = reference(lams.pop() if lams else 0.2)

    plot_separation(cav, bare, bend1, rabi1, out)
    plot_spectra(cav, bare, out)
    with open(out / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS)
        for r in bare + cav:
            w.writerow([r["run"], r["R"], r["lam"], r["steps"], f"{r['res']:.1f}", f"{r['bend']:.0f}",
                        f"{r['lp']:.0f}" if r["cavity"] else "", f"{r['up']:.0f}" if r["cavity"] else "",
                        f"{r['up'] - r['lp']:.0f}" if r["resolved"] else "", r["resolved"],
                        f"{r['rot']:.1f}", f"{r['drift']:.3f}", f"{r['sep_min']:.2f}"])

    print(f"single molecule: bend {bend1:.0f} cm^-1, Rabi {rabi1:.0f} cm^-1 "
          f"-> sqrt(2) x = {np.sqrt(2) * rabi1:.0f} cm^-1")
    print(f"{'run':40s} {'R':>5s} {'bend':>6s} {'LP':>6s} {'UP':>6s} {'Rabi':>6s} {'ratio':>6s} "
          f"{'rot':>6s} {'drift':>6s} {'sepmin':>6s}")
    for r in bare + cav:
        rabi = r["up"] - r["lp"] if r["resolved"] else np.nan
        print(f"{r['run']:40s} {r['R']:5.1f} {r['bend']:6.0f} "
              f"{(r['lp'] if r['cavity'] else np.nan):6.0f} {(r['up'] if r['cavity'] else np.nan):6.0f} "
              f"{rabi:6.0f} {rabi / rabi1:6.2f} {r['rot']:6.1f} {r['drift']:6.3f} {r['sep_min']:6.2f}")
    print(f"\nwrote {out}/separation.png, ir_spectra.png, peaks.csv")


if __name__ == "__main__":
    main()
