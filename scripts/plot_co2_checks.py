"""
Plot the CO2 single-molecule checks (scripts/campaigns/co2_checks.sh), NEP chi off,
lambda = 0.1, 10000 steps (6.7 cm^-1 resolution):

  orientation : Rabi splitting vs angle theta between cavity polarization and the
                CO2 axis (x), vs Omega_R(0) cos(theta). The molecule's maximum
                rotation during the run is read from md.traj (local only) and
                cached in peaks.csv so a fresh clone can still plot it.
  detuning    : LP / UP vs cavity frequency, with the exact linear model
                (one vibration + one cavity mode + dipole self-energy), one fit
                parameter Lambda = lam^2 (dmu/dQ)^2:
                  w±^2 = [(wv^2 + L + wc^2) ± sqrt((wv^2 + L - wc^2)^2 + 4 L wc^2)] / 2
  amplitude   : Rabi splitting vs initial C displacement (linear regime check).

Peaks: same method as plot_co2_sweep.py (the group benchmark's).

Usage:
    python scripts/plot_co2_checks.py
Writes to results/co2/2026-10-02_checks/:
    orientation.png, detuning.png, detuning_spectra.png, amplitude.png, peaks.csv
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_co2_sweep as sw  # noqa: E402  (spectrum, peaks, styling)

REPO = sw.REPO
RUNS = REPO / "runs" / "co2"
OUT = REPO / "results" / "co2" / "2026-10-02_checks"
BARE = RUNS / "2026-10-02_lambda_sweep" / "nep_bare"
CM_PER_AU = 219474.63
STYLE = sw.SERIES["nep_nochi"]          # same entity, same color/marker as the sweep
INK, INK2, GRID, SURFACE = sw.INK, sw.INK2, sw.GRID, sw.SURFACE


def dipole_file(d):
    return d / "dipole.dat" if (d / "dipole.dat").exists() else d / "dipole.dat.gz"


def max_rotation(d):
    """Largest angle (deg) the O-O axis turns away from its start, or None without md.traj."""
    if not (d / "md.traj").exists():
        return None
    from ase.io import read
    traj = read(d / "md.traj", "::200")
    ax = [t.positions[0] - t.positions[2] for t in traj]
    ax = [v / np.linalg.norm(v) for v in ax]
    return max(np.degrees(np.arccos(np.clip(abs(v @ ax[0]), -1, 1))) for v in ax)


def load(d):
    cfg = json.loads((d / "in.json").read_text())
    dip = np.loadtxt(dipole_file(d))
    freq, amp = sw.spectrum(dip[:, 2], dip[1, 1] - dip[0, 1])
    lp, up = sw.polariton_peaks(freq, amp)
    vec = np.asarray(cfg.get("lambda_vector", [[1, 0, 0]])[0], dtype=float)
    geom = REPO / "geometries" / Path(cfg["xyz_file"]).name     # in.json paths are absolute
    disp = float(np.loadtxt(geom, skiprows=2, usecols=1)[1])     # C x-coordinate
    return dict(name=d.name, freq=freq, amp=amp, lp=lp, up=up, rabi=up - lp,
                theta=np.degrees(np.arccos(abs(vec[0]) / np.linalg.norm(vec))),
                omega_c=cfg["omega_photon"][0] * CM_PER_AU if cfg.get("photons") else None, disp=disp)


def load_campaign(name, key):
    runs = sorted((load(d) for d in (RUNS / name).iterdir() if (d / "in.json").exists()), key=lambda r: r[key])
    return runs


def hopfield(L, wc, wv):
    a = wv ** 2 + L + wc ** 2
    b = np.sqrt((wv ** 2 + L - wc ** 2) ** 2 + 4 * L * wc ** 2)
    return np.sqrt((a - b) / 2), np.sqrt((a + b) / 2)


def points(ax, x, y, label=None):
    ax.plot(x, y, ls="none", marker=STYLE["marker"], ms=8, color=STYLE["color"], mec=SURFACE, mew=2,
            label=label, zorder=3)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cache = {}
    if (OUT / "peaks.csv").exists():
        cache = {r["run"]: r for r in csv.DictReader(open(OUT / "peaks.csv"))}
    bare = load(BARE)
    wv = bare["lp"]
    rows = []

    # ── orientation ──
    orient = load_campaign("2026-10-02_orientation", "theta")
    for r in orient:
        rot = max_rotation(RUNS / "2026-10-02_orientation" / r["name"])
        if rot is None and cache.get(r["name"], {}).get("max_rotation_deg"):
            rot = float(cache[r["name"]]["max_rotation_deg"])
        r["rot"] = rot
    r0 = orient[0]["rabi"]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    th = np.linspace(0, 90, 200)
    ax.plot(th, r0 * np.cos(np.radians(th)), color=INK2, lw=1, ls="--",
            label=rf"$\Omega_R(0°)\,\cos\theta$  ({r0:.0f} cm⁻¹ at 0°)", zorder=1)
    points(ax, [r["theta"] for r in orient], [r["rabi"] for r in orient], "NEP, χ neglected, λ = 0.1")
    for r in orient:
        if r["rot"] and r["rot"] > 1:
            ax.annotate(f"rotates {r['rot']:.0f}°", (r["theta"], r["rabi"]), xytext=(6, 8),
                        textcoords="offset points", fontsize=8, color=INK2)
    ax.set_xlabel("Angle θ between cavity polarization and CO₂ axis (°)")
    ax.set_ylabel(r"Rabi splitting $\Omega_R$ (cm⁻¹)")
    ax.set_title("CO₂ orientation check: follows cos θ to 45°, then the molecule turns\n"
                 "away from the polarization and the splitting collapses", loc="left", fontsize=10)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_xlim(-3, 93)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "orientation.png", dpi=160)
    plt.close(fig)

    # ── detuning ──
    det = load_campaign("2026-10-02_detuning", "omega_c")
    wc = np.array([r["omega_c"] for r in det])
    lp, up = np.array([r["lp"] for r in det]), np.array([r["up"] for r in det])
    fit = least_squares(lambda p: np.r_[hopfield(p[0], wc, wv)[0] - lp, hopfield(p[0], wc, wv)[1] - up], [1e5])
    L = fit.x[0]
    m_lp, m_up = hopfield(L, wc, wv)
    rms = np.sqrt(np.mean(np.r_[m_lp - lp, m_up - up] ** 2))
    x = np.linspace(2050, 2850, 300)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(x, x, color=INK2, lw=0.8, ls=":", label=r"uncoupled cavity ($\omega = \omega_c$)")
    ax.axhline(wv, color=INK2, lw=0.8, ls="-.", label=f"bare stretch ({wv:.0f} cm⁻¹)")
    ax.axhline(np.sqrt(wv ** 2 + L), color=INK2, lw=0.8, ls="--",
               label=rf"self-energy-shifted stretch $\sqrt{{\omega_v^2+\Lambda}}$ = {np.sqrt(wv ** 2 + L):.0f}")
    for i, y in enumerate(hopfield(L, x, wv)):
        ax.plot(x, y, color=INK, lw=1.2, zorder=1,
                label=rf"exact linear model, fit $\sqrt{{\Lambda}}$ = {np.sqrt(L):.0f} cm⁻¹ (rms {rms:.1f})" if i == 0 else None)
    points(ax, wc, lp, "NEP, χ neglected, λ = 0.1")
    points(ax, wc, up)
    ax.set_xlabel(r"Cavity frequency $\omega_c$ (cm⁻¹)")
    ax.set_ylabel("Polariton frequency (cm⁻¹)")
    ax.set_title("CO₂ detuning scan: anti-crossing matches the exact linear model", loc="left", fontsize=11)
    ax.set_xlim(2050, 2850)
    ax.set_ylim(1950, 2950)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "detuning.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    for i, r in enumerate(det):
        f, a = r["freq"], r["amp"]
        m = (f > 1900) & (f < 3000)
        ax.plot(f[m], a[m] / a[m].max() * 0.85 + i, color=STYLE["color"], lw=1.3)
        ax.text(1910, i + 0.55, rf"$\omega_c$ = {r['omega_c']:.0f}", fontsize=9)
        ax.plot([r["omega_c"]] * 2, [i, i + 0.92], color=INK2, lw=1, ls="--")
    ax.axvline(wv, color=INK2, lw=0.8, ls="-.")
    ax.set_xlim(1900, 3000)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Frequency (cm⁻¹)")
    ax.set_title("CO₂ spectra vs cavity frequency (dashed = ω_c, dash-dot = bare stretch)", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "detuning_spectra.png", dpi=160)
    plt.close(fig)

    # ── amplitude ──
    amp = load_campaign("2026-10-02_amplitude", "disp")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    points(ax, [r["disp"] for r in amp], [r["rabi"] for r in amp], "NEP, χ neglected, λ = 0.1")
    ax.set_xscale("log")
    ax.set_xticks([r["disp"] for r in amp])
    ax.set_xticklabels([f"{r['disp']:g}" for r in amp])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Initial C displacement along the axis (Å)")
    ax.set_ylabel(r"Rabi splitting $\Omega_R$ (cm⁻¹)")
    ax.set_ylim(0, max(r["rabi"] for r in amp) * 1.3)
    ax.set_title("CO₂ amplitude check: splitting independent of starting amplitude", loc="left", fontsize=10)
    ax.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "amplitude.png", dpi=160)
    plt.close(fig)

    with open(OUT / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["check", "run", "theta_deg", "omega_c_cm-1", "C_disp_A", "LP_cm-1", "UP_cm-1", "Rabi_cm-1", "max_rotation_deg"])
        for check, runs in (("orientation", orient), ("detuning", det), ("amplitude", amp)):
            for r in runs:
                rot = r.get("rot")
                w.writerow([check, r["name"], f"{r['theta']:.0f}", f"{r['omega_c']:.0f}", r["disp"],
                            f"{r['lp']:.0f}", f"{r['up']:.0f}", f"{r['rabi']:.0f}", "" if rot is None else f"{rot:.1f}"])
    print(f"detuning fit: sqrt(Lambda) = {np.sqrt(L):.0f} cm-1, rms {rms:.1f} cm-1; wrote {OUT.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
