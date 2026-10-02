"""
Plot the H2O single-molecule checks (scripts/campaigns/h2o_checks.sh):

  orientation : Rabi splitting vs angle theta between the cavity polarization and
                the bend dipole (x), compared with Omega_R(0) * cos(theta).
                CAVEAT (found 2026-10-02): for 0 < theta < 90 the cavity field
                exerts a torque on H2O's permanent dipole and the molecule rotates
                (49 deg at theta=30, 113 deg at theta=60 within 400 fs), so this is
                not a clean cos(theta) test for a free polar molecule.
  detuning    : LP / UP vs cavity frequency, compared with the two-coupled-
                oscillator anti-crossing  w± = (wc + wv)/2 ± sqrt(g^2 + (wc - wv)^2/4),
                wv = bare bend (from the lambda-sweep bare run), g = Omega_R/2 on
                resonance (from the orientation theta = 0 run, same lambda).
                The data show the anti-crossing but sit ~40-60 cm^-1 below this
                guide (and below the exact linear model with dipole self-energy),
                so the curve is a guide, not a fit.

Spectra and peak finding come from plot_h2o_sweep.py (unwindowed, zero-padded FFT
of dipole_x). For the detuning scan the far-detuned polariton is weak, so the
second peak only needs 10% of the first, but must sit >= 2 resolution bins away
(the finite-length side lobes sit closer than that).

Usage:
    python scripts/plot_h2o_checks.py
Writes to results/h2o/2026-10-02_checks/:
    orientation.png, detuning.png, orientation_spectra.png, detuning_spectra.png, peaks.csv
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_h2o_sweep as sw  # noqa: E402  (spectrum, peaks, styling)

REPO = sw.REPO
RUNS = REPO / "runs" / "h2o"
ORIENT, DETUNE = RUNS / "2026-10-02_orientation", RUNS / "2026-10-02_detuning"
BARE = RUNS / "2026-10-02_lambda_sweep" / "pyscf_bare"
OUT = REPO / "results" / "h2o" / "2026-10-02_checks"
CM_PER_AU = 219474.63


def load(run):
    cfg = json.loads((run / "in.json").read_text())
    f = run / "dipole.dat" if (run / "dipole.dat").exists() else run / "dipole.dat.gz"
    dip = np.loadtxt(f)
    freq, amp = sw.spectrum(dip[:, 2], cfg["timestep"])
    res = 1 / (cfg["steps"] * cfg["timestep"] * 1e-15) / sw.C_CM_S
    vec = np.asarray(cfg.get("lambda_vector", [[1, 0, 0]])[0], dtype=float)
    return dict(name=run.name, cfg=cfg, freq=freq, amp=amp, res=res,
                theta=np.degrees(np.arccos(abs(vec[0]) / np.linalg.norm(vec))),
                omega_c=cfg["omega_photon"][0] * CM_PER_AU if cfg.get("photons") else None)


def two_peaks(r, ratio, min_sep_bins=0.0):
    """(LP, UP) or None if not split by this rule."""
    peaks = sw.local_maxima(r["freq"], r["amp"], sw.BEND_BAND)
    first = peaks[0]
    for f, a in peaks[1:]:
        if a >= ratio * first[1] and abs(f - first[0]) >= min_sep_bins * r["res"]:
            return tuple(sorted((first[0], f))) + (min(a, first[1]) / max(a, first[1]),)
    return None


def stacked(ax, runs, label_of, omega_of=None):
    for i, r in enumerate(runs):
        f, a = r["freq"], r["amp"]
        m = (f > sw.BEND_BAND[0]) & (f < sw.BEND_BAND[1])
        ax.plot(f[m], a[m] / a[m].max() * 0.85 + i, color=sw.COLOR, lw=1.5)
        ax.text(sw.BEND_BAND[0] + 15, i + 0.55, label_of(r), fontsize=9, color=sw.INK)
        if r.get("lp"):
            for pk in (r["lp"], r["up"]):
                ax.plot([pk, pk], [i, i + 0.92], color=sw.INK2, lw=0.8, ls=":")
        if omega_of:
            ax.plot([omega_of(r)] * 2, [i, i + 0.92], color=sw.INK2, lw=1, ls="--")
    ax.set_xlim(*sw.BEND_BAND)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Frequency (cm⁻¹)")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bare = load(BARE)
    wv = sw.local_maxima(bare["freq"], bare["amp"], sw.BEND_BAND)[0][0]
    lam = None
    rows = []

    # ── orientation ──
    orient = sorted((load(d) for d in ORIENT.iterdir() if (d / "in.json").exists()), key=lambda r: r["theta"])
    for r in orient:
        pk = two_peaks(r, sw.SPLIT_RATIO)
        r["lp"], r["up"] = (pk[0], pk[1]) if pk else (None, None)
        lam = r["cfg"]["lambda_photon"][0]
        rows.append(["orientation", r["name"], f"{r['theta']:.0f}", f"{r['omega_c']:.0f}",
                     f"{r['lp']:.0f}" if pk else "", f"{r['up']:.0f}" if pk else "",
                     f"{r['up'] - r['lp']:.0f}" if pk else "unresolved"])
    rabi0 = orient[0]["up"] - orient[0]["lp"]

    fig, ax = plt.subplots(figsize=(7.5, 5))
    th = np.linspace(0, 90, 200)
    ax.plot(th, rabi0 * np.cos(np.radians(th)), color=sw.INK2, lw=1, ls="--",
            label=rf"$\Omega_R(0°)\,\cos\theta$  ({rabi0:.0f} cm⁻¹ at 0°)", zorder=1)
    ax.axhspan(0, orient[0]["res"], color=sw.GRID, alpha=0.6, lw=0, zorder=0)
    ax.text(88, orient[0]["res"] / 2, f"below resolution ({orient[0]['res']:.0f} cm⁻¹)",
            fontsize=8.5, color=sw.INK2, ha="right", va="center")
    ok = [r for r in orient if r["lp"]]
    ax.plot([r["theta"] for r in ok], [r["up"] - r["lp"] for r in ok], ls="none", marker=sw.MARKER,
            ms=9, color=sw.COLOR, mec=sw.SURFACE, mew=2, label=sw.LABEL, zorder=3)
    nores = [r for r in orient if not r["lp"]]
    if nores:
        ax.plot([r["theta"] for r in nores], [0] * len(nores), ls="none", marker=sw.MARKER, ms=9,
                mfc=sw.SURFACE, mec=sw.COLOR, mew=1.5, label="unresolved (single peak)", zorder=3)
    ax.set_xlabel("Angle θ between cavity polarization and bend dipole (°)")
    ax.set_ylabel(r"Rabi splitting $\Omega_R$ (cm⁻¹)")
    ax.set_title(f"H₂O orientation check (λ = {lam:g}): NOT a clean cos θ test.\n"
                 "At 0° < θ < 90° the cavity field torques the permanent dipole and the molecule rotates",
                 loc="left", fontsize=10)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_xlim(-3, 93)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "orientation.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6.5))
    stacked(ax, orient, lambda r: f"θ = {r['theta']:.0f}°" + ("" if r["lp"] else "  (unresolved)"))
    ax.axvline(wv, color=sw.INK2, lw=1, ls="--")
    ax.set_title(f"H₂O spectra vs cavity polarization angle (λ = {lam:g})", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "orientation_spectra.png", dpi=160)
    plt.close(fig)

    # ── detuning ──
    detune = sorted((load(d) for d in DETUNE.iterdir() if (d / "in.json").exists()), key=lambda r: r["omega_c"])
    for r in detune:
        pk = two_peaks(r, 0.10, min_sep_bins=2.0)
        r["lp"], r["up"] = (pk[0], pk[1]) if pk else (None, None)
        rows.append(["detuning", r["name"], "0", f"{r['omega_c']:.0f}",
                     f"{r['lp']:.0f}" if pk else "", f"{r['up']:.0f}" if pk else "",
                     f"{r['up'] - r['lp']:.0f}" if pk else "unresolved"])
    g = rabi0 / 2
    wc = np.linspace(1250, 1900, 300)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(wc, wc, color=sw.INK2, lw=0.8, ls=":", label="uncoupled cavity (ω = ω_c)")
    ax.axhline(wv, color=sw.INK2, lw=0.8, ls="-.", label=f"uncoupled bend ({wv:.0f} cm⁻¹)")
    for sgn in (-1, 1):
        ax.plot(wc, (wc + wv) / 2 + sgn * np.sqrt(g ** 2 + (wc - wv) ** 2 / 4), color=sw.INK, lw=1.2,
                label=f"simple coupled-oscillator guide (g = {g:.0f} cm⁻¹; not a fit)" if sgn < 0 else None, zorder=1)
    ok = [r for r in detune if r["lp"]]
    for key in ("lp", "up"):
        ax.plot([r["omega_c"] for r in ok], [r[key] for r in ok], ls="none", marker=sw.MARKER, ms=9,
                color=sw.COLOR, mec=sw.SURFACE, mew=2, label=sw.LABEL if key == "lp" else None, zorder=3)
    ax.set_xlabel(r"Cavity frequency $\omega_c$ (cm⁻¹)")
    ax.set_ylabel("Polariton frequency (cm⁻¹)")
    ax.set_title(f"H₂O detuning scan: anti-crossing (λ = {lam:g})", loc="left", fontsize=11)
    ax.set_xlim(1250, 1900)
    ax.set_ylim(1150, 2050)
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "detuning.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    stacked(ax, detune, lambda r: rf"$\omega_c$ = {r['omega_c']:.0f}", omega_of=lambda r: r["omega_c"])
    ax.axvline(wv, color=sw.INK2, lw=0.8, ls="-.")
    ax.set_title(f"H₂O spectra vs cavity frequency (λ = {lam:g}; dashed = ω_c, dash-dot = bare bend)",
                 loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "detuning_spectra.png", dpi=160)
    plt.close(fig)

    with open(OUT / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["check", "run", "theta_deg", "omega_c_cm-1", "LP_cm-1", "UP_cm-1", "Rabi_cm-1"])
        w.writerows(rows)
    for row in rows:
        print("  ".join(row))
    print(f"\nbare bend {wv:.0f} cm⁻¹; wrote {OUT.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
