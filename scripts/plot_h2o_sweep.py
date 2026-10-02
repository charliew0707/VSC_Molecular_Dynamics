"""
Plot an H2O lambda-sweep campaign: IR spectra, polariton branches, Rabi splitting.

Reads every run folder in a campaign (made by scripts/campaigns/h2o_lambda_sweep.sh):
the bare run (*_bare) is lambda = 0, cavity runs are swept in lambda.

Spectra: FFT of the mean-subtracted dipole_x (the component `infrared -m ase` uses;
the H2O geometry has its C2 axis, and so the bend dipole, along x), with no window
and zero-padded 8x. Zero-padding only interpolates between frequency points; the
true resolution is 1 / (steps * dt), ~83 cm^-1 for 800 x 0.5 fs. Windowing (e.g.
Hanning) would blur the small H2O splittings further at this trajectory length.

Polaritons: the two strongest local maxima in the bend band (1000-2300 cm^-1),
counted as split only if the weaker is >= 50% of the stronger (finite-length
side lobes are ~20%). Otherwise the run is marked "unresolved".

Usage:
    python scripts/plot_h2o_sweep.py runs/h2o/2026-10-02_lambda_sweep
Writes to results/h2o/<campaign>/:
    ir_spectra.png          stacked IR spectra vs lambda (bend + 2500-4200 cm^-1 region)
    polariton_branches.png  LP / UP vs lambda
    rabi_splitting.png      UP - LP vs lambda
    peaks.csv               every run's LP, UP, Rabi splitting
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
C_CM_S = 2.99792458e10
PAD = 8
HIGH_SCALE = 3                  # magnification of the 2500-4200 cm^-1 panel
BEND_BAND = (1000.0, 2300.0)
HIGH_BAND = (2500.0, 4200.0)     # symmetric stretch + polariton overtones/combinations
SPLIT_RATIO = 0.5
# Earlier result: runs/h2o/pyscf_cavity (2000 steps, z-oriented geometry, Hanning analysis)
EARLIER = dict(lam=0.1, lp=1517.0, up=1642.0)

# Same entity, same color as the CO2 sweep: PySCF, chi neglected = categorical slot 3
COLOR, MARKER, LABEL = "#1baf7a", "^", "PySCF (LDA/cc-pVDZ), χ neglected, 800 steps"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "font.size": 10, "legend.frameon": False,
})


# ── data ──────────────────────────────────────────────────────────────────────

def spectrum(dipole_x, dt_fs):
    s = dipole_x - dipole_x.mean()
    amp = np.abs(np.fft.rfft(s, n=PAD * len(s)))
    freq = np.fft.rfftfreq(PAD * len(s), d=dt_fs * 1e-15) / C_CM_S
    return freq, amp


def local_maxima(freq, amp, band):
    m = (freq > band[0]) & (freq < band[1])
    f, a = freq[m], amp[m]
    k = [i for i in range(1, len(a) - 1) if a[i] > a[i - 1] and a[i] > a[i + 1]]
    return sorted(((f[i], a[i]) for i in k), key=lambda p: -p[1])


def polaritons(freq, amp):
    """(LP, UP, resolved). Unresolved: LP = UP = the single bend peak."""
    peaks = local_maxima(freq, amp, BEND_BAND)
    if len(peaks) >= 2 and peaks[1][1] >= SPLIT_RATIO * peaks[0][1]:
        lo, hi = sorted((peaks[0][0], peaks[1][0]))
        return lo, hi, True
    return peaks[0][0], peaks[0][0], False


def load_campaign(campaign):
    def dipole_file(d):
        return d / "dipole.dat" if (d / "dipole.dat").exists() else d / "dipole.dat.gz"

    runs = []
    for d in sorted(p for p in campaign.iterdir() if dipole_file(p).exists()):
        cfg = json.loads((d / "in.json").read_text())
        dip = np.loadtxt(dipole_file(d))
        if len(dip) < cfg["steps"] + 1:
            print(f"skipping {d.name}: incomplete ({len(dip)}/{cfg['steps'] + 1} steps)")
            continue
        lam = cfg["lambda_photon"][0] if cfg.get("photons") else 0.0
        freq, amp = spectrum(dip[:, 2], cfg["timestep"])
        lp, up, resolved = polaritons(freq, amp)
        runs.append(dict(name=d.name, lam=lam, steps=cfg["steps"], freq=freq, amp=amp,
                         lp=lp, up=up, resolved=resolved or lam == 0,
                         res=1 / (cfg["steps"] * cfg["timestep"] * 1e-15) / C_CM_S,
                         omega_c=cfg["omega_photon"][0] * 219474.63 if cfg.get("photons") else None))
    return sorted(runs, key=lambda r: r["lam"])


# ── figures ───────────────────────────────────────────────────────────────────

def plot_spectra(runs, out):
    omega_c = next(r["omega_c"] for r in runs if r["omega_c"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 7), sharey=True,
                                   gridspec_kw=dict(width_ratios=[1.3, 1]))
    for i, r in enumerate(runs):
        color = COLOR if r["lam"] > 0 else INK2
        f, a = r["freq"], r["amp"]
        bend_max = a[(f > BEND_BAND[0]) & (f < BEND_BAND[1])].max()
        for ax, band, scale in ((ax1, BEND_BAND, 1), (ax2, HIGH_BAND, HIGH_SCALE)):
            m = (f > band[0]) & (f < band[1])
            ax.plot(f[m], a[m] / bend_max * 0.85 * scale + i, color=color, lw=1.5)
        label = "no cavity" if r["lam"] == 0 else f"λ = {r['lam']:g}"
        if r["lam"] > 0 and not r["resolved"]:
            label += "  (unresolved)"
        ax1.text(BEND_BAND[0] + 15, i + 0.55, label, fontsize=9, color=INK)
        if r["lam"] > 0 and r["resolved"]:
            for pk in (r["lp"], r["up"]):
                ax1.plot([pk, pk], [i, i + 0.92], color=INK2, lw=0.8, ls=":")
    ax1.axvline(omega_c, color=INK2, lw=1, ls="--")
    ax1.text(omega_c + 15, 1.0, r"$\omega_c$", color=INK2, fontsize=10,
             transform=ax1.get_xaxis_transform(), va="top")
    ax1.set_title("Bend region (dotted = LP / UP)", loc="left", fontsize=11)
    ax2.set_title(f"Sym. stretch + polariton overtones (×{HIGH_SCALE})", loc="left", fontsize=11)
    for ax, band in ((ax1, BEND_BAND), (ax2, HIGH_BAND)):
        ax.set_xlim(*band)
        ax.set_xlabel("Frequency (cm⁻¹)")
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
    ax2.text(3840, 0.95, "sym. stretch\n(far from $\\omega_c$,\nnot split)", fontsize=8.5, color=INK2, va="top")
    ax2.text(2520, len(runs) - 0.25, "LP+UP, 2LP, 2UP (anharmonic overtones)\ngrow with λ", fontsize=8.5,
             color=INK2, va="top")
    ax1.set_ylabel("IR intensity (normalized to the bend, offset by λ)")
    fig.suptitle(rf"H₂O IR spectra in a cavity resonant with the bend ($\omega_c$ = {omega_c:.0f} cm⁻¹), "
                 f"PySCF, {runs[0]['steps']} steps ({runs[0]['res']:.0f} cm⁻¹ resolution)",
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out / "ir_spectra.png", dpi=160)
    plt.close(fig)


def _earlier_point(ax, value, label):
    ax.plot(EARLIER["lam"], value, ls="none", marker="D", ms=13, mfc="none", mec=INK, mew=1.3,
            label=label, zorder=4)


def plot_branches(runs, out):
    omega_c = next(r["omega_c"] for r in runs if r["omega_c"])
    bare = next(r for r in runs if r["lam"] == 0)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    res = [r for r in runs if r["resolved"] and r["lam"] > 0]
    unres = [r for r in runs if not r["resolved"]]
    for key in ("lp", "up"):
        ax.plot([0] + [r["lam"] for r in res], [bare["lp"]] + [r[key] for r in res],
                color=COLOR, lw=1.2, alpha=0.6, zorder=2)
        ax.plot([r["lam"] for r in res], [r[key] for r in res], ls="none", marker=MARKER, ms=9,
                color=COLOR, mec=SURFACE, mew=2, label=LABEL if key == "lp" else None, zorder=3)
    ax.plot([0], [bare["lp"]], ls="none", marker="o", ms=8, color=INK2, mec=SURFACE, mew=2,
            label=f"no cavity: bend at {bare['lp']:.0f} cm⁻¹", zorder=3)
    if unres:
        ax.plot([r["lam"] for r in unres], [r["lp"] for r in unres], ls="none", marker=MARKER, ms=9,
                mfc=SURFACE, mec=COLOR, mew=1.5, label="unresolved (splitting < resolution)", zorder=3)
    _earlier_point(ax, EARLIER["lp"], "earlier run: 2000 steps, Hanning (runs/h2o/pyscf_cavity)")
    _earlier_point(ax, EARLIER["up"], None)
    ax.axhline(omega_c, color=INK2, lw=0.8, ls=":")
    ax.text(0.316, omega_c, r"$\omega_c$", color=INK2, fontsize=10, va="center")
    ax.text(0.305, res[-1]["up"], "UP", color=INK2, fontsize=10, va="center")
    ax.text(0.305, res[-1]["lp"], "LP", color=INK2, fontsize=10, va="center")
    ax.legend(fontsize=8.5, loc="lower left")
    ax.set_xlabel("Coupling strength λ (a.u.)")
    ax.set_ylabel("Polariton frequency (cm⁻¹)")
    ax.set_title("H₂O bend polariton branches vs coupling", loc="left", fontsize=11)
    ax.set_xlim(-0.005, 0.33)
    fig.tight_layout()
    fig.savefig(out / "polariton_branches.png", dpi=160)
    plt.close(fig)


def plot_rabi(runs, out):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    res = [r for r in runs if r["resolved"] and r["lam"] > 0]
    lam = np.array([r["lam"] for r in res])
    rabi = np.array([r["up"] - r["lp"] for r in res])
    slope = (lam @ rabi) / (lam @ lam)                     # least-squares line through 0
    x = np.linspace(0, 0.32, 50)
    ax.plot(x, slope * x, color=INK2, lw=1, ls="--", label=rf"linear fit through 0: $\Omega_R$ ≈ {slope:.0f}·λ cm⁻¹", zorder=1)
    ax.axhspan(0, runs[0]["res"], color=GRID, alpha=0.6, lw=0, zorder=0)
    ax.text(0.31, runs[0]["res"] / 2, f"below resolution ({runs[0]['res']:.0f} cm⁻¹)",
            fontsize=8.5, color=INK2, ha="right", va="center")
    ax.plot(lam, rabi, ls="none", marker=MARKER, ms=9, color=COLOR, mec=SURFACE, mew=2, label=LABEL, zorder=3)
    _earlier_point(ax, EARLIER["up"] - EARLIER["lp"], "earlier run: 2000 steps, Hanning (runs/h2o/pyscf_cavity)")
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_xlabel("Coupling strength λ (a.u.)")
    ax.set_ylabel(r"Rabi splitting $\Omega_R$ = UP − LP (cm⁻¹)")
    ax.set_title("H₂O bend Rabi splitting vs coupling", loc="left", fontsize=11)
    ax.set_xlim(-0.005, 0.32)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out / "rabi_splitting.png", dpi=160)
    plt.close(fig)


def write_table(runs, out):
    with open(out / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "lambda", "steps", "resolution_cm-1", "resolved", "LP_cm-1", "UP_cm-1", "Rabi_cm-1"])
        for r in runs:
            w.writerow([r["name"], r["lam"], r["steps"], f"{r['res']:.1f}", r["resolved"],
                        f"{r['lp']:.0f}", f"{r['up']:.0f}", f"{r['up'] - r['lp']:.0f}" if r["resolved"] else ""])


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    campaign = Path(sys.argv[1]).resolve()
    out = REPO / "results" / "h2o" / campaign.name
    out.mkdir(parents=True, exist_ok=True)
    runs = load_campaign(campaign)
    plot_spectra(runs, out)
    plot_branches(runs, out)
    plot_rabi(runs, out)
    write_table(runs, out)
    print(f"{LABEL}  ({runs[0]['res']:.0f} cm⁻¹ resolution)")
    for r in runs:
        tail = f"Ω_R {r['up'] - r['lp']:5.0f} cm⁻¹" if r["resolved"] and r["lam"] else ("bare bend" if not r["lam"] else "unresolved")
        print(f"  λ = {r['lam']:<5g}  LP {r['lp']:6.0f}  UP {r['up']:6.0f}  {tail}")
    print(f"\nwrote {out.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
