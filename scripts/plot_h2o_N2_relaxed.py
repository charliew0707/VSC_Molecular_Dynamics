"""
Figures and peak table for the 2026-11-18 two-water work (relaxed photon start,
head-to-tail vs side-by-side, long runs). Used by reports/2026-11-18_weekly_report.tex.

Writes to results/h2o_N2/2026-11-18_relaxed_photon/:
    arrangements.png   the two starting geometries seen from above (x-y plane)
    photon_start.png   q = 0 vs relaxed photon start: x-dipole and splitting
    spectra.png        bend region: 1 water, 2 waters (both arrangements), bare pair
    motion.png         separation and libration angle vs time (long runs)
    ratio.png          Rabi(2 waters) / Rabi(1 water) for every comparison, vs sqrt(2)
    peaks.csv          LP / UP / Rabi for every run and window used
    motion.csv         separation and axis angle vs time (cached: position.dat stays local)

Peaks: zero-padded FFT of the x dipole and the >= 50% height rule from
plot_h2o_sweep.py. Head-to-tail long run: first 1 ps only (the waters collide later),
compared with the single water over the same 1 ps.

Usage:
    python scripts/plot_h2o_N2_relaxed.py
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import io as aio

import plot_h2o_sweep as sw

REPO = sw.REPO
OUT = REPO / "results" / "h2o_N2" / "2026-11-18_relaxed_photon"
INK, INK2 = sw.INK, sw.INK2
GREEN, BLUE, ORANGE, GOLD = "#1baf7a", "#2a78d6", "#eb6834", "#eda100"
GREY = "#9a9893"
AU_FS = 0.0241888432658

R = lambda p: REPO / "runs" / p
RUNS = {
    "1w_q0":        R("h2o/2026-10-02_lambda_sweep/pyscf_cav_w1584_lam0.2_nochi"),
    "1w_rel":       R("h2o/2026-11-18_relaxed_photon/pyscf_cav_w1584_lam0.2_relaxed"),
    "1w_rel_long":  R("h2o/2026-11-18_relaxed_photon_long/pyscf_cav_w1584_lam0.2_relaxed"),
    "ht_q0":        R("h2o_N2/2026-11-18_separation_scan/pyscf_cav_w1584_lam0.2_R8.0"),
    "ht_bare":      R("h2o_N2/2026-11-18_separation_scan/pyscf_R8.0_bare"),
    "ht_rel":       R("h2o_N2/2026-11-18_relaxed_photon/pyscf_cav_w1584_lam0.2_R8.0_relaxed"),
    "ht_rel_long":  R("h2o_N2/2026-11-18_relaxed_photon_long/pyscf_cav_w1584_lam0.2_R8.0_relaxed"),
    "sbs_rel":      R("h2o_N2/2026-11-18_sidebyside/pyscf_cav_w1584_lam0.2_sidebyside_R8.0_relaxed"),
    "sbs_bare":     R("h2o_N2/2026-11-18_sidebyside/pyscf_sidebyside_R8.0_bare"),
    "sbs_rel_long": R("h2o_N2/2026-11-18_sidebyside_long/pyscf_cav_w1584_lam0.2_sidebyside_R8.0_relaxed"),
}
FIRST_PS = 2000      # steps in the first 1 ps (dt = 0.5 fs)


def dipole(run):
    f = run / "dipole.dat"
    d = np.loadtxt(f if f.exists() else run / "dipole.dat.gz")
    return (d[:, 1] - d[0, 1]) * AU_FS, d[:, 2]          # time (fs), mu_x (a.u.)


def peaks(run, stop=None):
    cfg = json.loads((run / "in.json").read_text())
    t, mu = dipole(run)
    mu = mu[:stop] if stop else mu
    freq, amp = sw.spectrum(mu, cfg["timestep"])
    lp, up, resolved = sw.polaritons(freq, amp)
    steps = len(mu) - 1
    return dict(lp=lp, up=up, rabi=up - lp if resolved else None, steps=steps,
                res=1 / (steps * cfg["timestep"] * 1e-15) / sw.C_CM_S, freq=freq, amp=amp)


# ── motion (position.dat on this machine, else the cached csv) ───────────────

def axes_and_centres(run):
    cfg = json.loads((run / "in.json").read_text())
    start = aio.read(REPO / "geometries" / Path(cfg["xyz_file"]).name)
    n = len(start)
    d = np.loadtxt(run / "position.dat")[::10]
    pos = d[:, 2:2 + 3 * n].reshape(len(d), n, 3)
    m = start.get_masses()
    out = []
    for sl in (slice(0, 3), slice(3, 6)):
        p = pos[:, sl]
        ax = p[:, 1:].mean(1) - p[:, 0]
        ax /= np.linalg.norm(ax, axis=1)[:, None]
        out.append((ax, (p * m[sl][None, :, None]).sum(1) / m[sl].sum()))
    t = d[:, 0] * 0.5                                      # step -> fs
    sep = np.linalg.norm(out[0][1] - out[1][1], axis=1)
    # signed in-plane angle of each C2 axis from its starting direction (-x)
    ang = [np.degrees(np.arctan2(a[:, 1], -a[:, 0])) for a, _ in out]
    return t, sep, ang[0], ang[1]


def motion():
    cache = OUT / "motion.csv"
    keys = ("ht_rel_long", "sbs_rel_long")
    if all((RUNS[k] / "position.dat").exists() for k in keys):
        rows = []
        for k in keys:
            t, sep, a1, a2 = axes_and_centres(RUNS[k])
            rows += [dict(run=k, time_fs=f"{a:.1f}", separation_A=f"{b:.4f}",
                          angle_mol1_deg=f"{c:.3f}", angle_mol2_deg=f"{d:.3f}")
                     for a, b, c, d in zip(t, sep, a1, a2)]
        with open(cache, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    data = {}
    with open(cache) as f:
        for r in csv.DictReader(f):
            data.setdefault(r["run"], []).append([float(r["time_fs"]), float(r["separation_A"]),
                                                  float(r["angle_mol1_deg"]), float(r["angle_mol2_deg"])])
    return {k: np.array(v) for k, v in data.items()}


# ── figures ──────────────────────────────────────────────────────────────────

def draw_water(ax, atoms, colour):
    p = atoms.positions
    for i in range(0, len(p), 3):
        o, h1, h2 = p[i], p[i + 1], p[i + 2]
        for h in (h1, h2):
            ax.plot([o[0], h[0]], [o[1], h[1]], color=INK2, lw=2, zorder=1)
        ax.scatter([o[0]], [o[1]], s=180, color="#d64545", zorder=2, edgecolor=INK, lw=0.6)
        ax.scatter([h1[0], h2[0]], [h1[1], h2[1]], s=70, color="white", zorder=2, edgecolor=INK, lw=0.6)
        d = (h1 + h2) / 2 - o
        d /= np.linalg.norm(d)
        base = o[:2] + np.array([0, -1.1])
        ax.annotate("", xy=base + 1.6 * d[:2], xytext=base,
                    arrowprops=dict(arrowstyle="-|>", color=colour, lw=2))


def plot_arrangements():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    for ax, name, colour, title, note in [
        (axs[0], "h2o_N2_headtail-x_R8.0A.xyz", BLUE, "Head to tail (shifted along x)",
         "dipoles in a line: attract"),
        (axs[1], "h2o_N2_sidebyside-y_R8.0A.xyz", ORANGE, "Side by side (shifted along y)",
         "dipoles parallel, side by side: repel")]:
        atoms = aio.read(REPO / "geometries" / name)
        atoms.positions -= atoms.positions.mean(0)
        draw_water(ax, atoms, colour)
        ax.set_aspect("equal")
        ax.set_xlim(-7, 7)
        ax.set_ylim(-6.5, 6.5)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_title(title, fontsize=10)
        ax.text(0, -5.8, note, ha="center", fontsize=9, color=INK2)
        ax.annotate("", xy=(6.3, 5.6), xytext=(3.8, 5.6),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.5))
        ax.text(5.05, 4.6, "cavity field", ha="center", fontsize=8, color=INK)
    fig.text(0.5, 0.005, "R = 8 Å between molecule centres; coloured arrows = each molecule's dipole (along the C$_2$ axis, O towards the H atoms)",
             ha="center", fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "arrangements.png", dpi=200)
    plt.close(fig)


def plot_photon_start(pk):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw=dict(width_ratios=[1.4, 1]))
    ax = axs[0]
    for k, c, lab in [("1w_q0", GREY, "photon starts at q = 0"), ("1w_rel", GREEN, "relaxed start, q = λμ/ω")]:
        t, mu = dipole(RUNS[k])
        ax.plot(t, (mu - mu[0]) * 1000, color=c, lw=1.4, label=lab)
    ax.set_xlim(0, 400)
    ax.set_xlabel("time (fs)")
    ax.set_ylabel("change in x dipole (10$^{-3}$ a.u.)")
    ax.set_title("One water, λ = 0.2: how hard the start shakes the molecule", fontsize=10)
    ax.legend(loc="lower right", fontsize=9)

    ax = axs[1]
    groups = [("one water", "1w_q0", "1w_rel"), ("two waters\n(head to tail)", "ht_q0", "ht_rel")]
    w = 0.36
    for i, (lab, a, b) in enumerate(groups):
        for j, (k, c, l) in enumerate([(a, GREY, "q = 0"), (b, GREEN, "relaxed")]):
            v = pk[k]["rabi"]
            ax.bar(i + (j - 0.5) * w, v, w, color=c, label=l if i == 0 else None)
            ax.text(i + (j - 0.5) * w, v + 6, f"{v:.0f}", ha="center", fontsize=9)
    ax.set_xticks([0, 1], [g[0] for g in groups])
    ax.set_ylabel("Rabi splitting (cm$^{-1}$)")
    ax.set_ylim(0, 440)
    ax.set_title("Splitting, 800 steps", fontsize=10)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(OUT / "photon_start.png", dpi=200)
    plt.close(fig)


def plot_spectra(pk):
    rows = [("1w_rel_long", GREEN, "1 water, cavity (2 ps)"),
            ("ht_rel_long_1ps", BLUE, "2 waters head to tail, cavity (first 1 ps)"),
            ("sbs_rel_long", ORANGE, "2 waters side by side, cavity (2 ps)"),
            ("sbs_bare", GREY, "2 waters side by side, no cavity (0.4 ps)")]
    fig, ax = plt.subplots(figsize=(8, 5.8))
    for i, (k, c, lab) in enumerate(rows):
        p = pk[k]
        m = (p["freq"] > 1100) & (p["freq"] < 2050)
        a = p["amp"][m] / p["amp"][m].max()
        off = (len(rows) - 1 - i) * 1.4
        ax.plot(p["freq"][m], a + off, color=c, lw=1.6)
        ax.text(1105, off + 1.1, lab, fontsize=9, color=c, va="bottom")
        if p["rabi"]:
            for f in (p["lp"], p["up"]):
                ax.plot([f], [off + 1.03], marker="v", color=c, ms=6)
            ax.text(p["up"] + 25, off + 0.45, f"Ω = {p['rabi']:.0f} cm$^{{-1}}$", fontsize=9, color=c)
    ax.axvline(1584, color=INK2, ls=":", lw=1)
    ax.text(1590, len(rows) * 1.4 - 0.1, "cavity = bare bend, 1584", fontsize=8, color=INK2)
    ax.set_xlim(1100, 2050)
    ax.set_yticks([])
    ax.set_xlabel("frequency (cm$^{-1}$)")
    ax.set_ylabel("IR intensity (each normalised, offset)")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(OUT / "spectra.png", dpi=200)
    plt.close(fig)


def plot_motion(mo):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    ht, sbs = mo["ht_rel_long"], mo["sbs_rel_long"]
    ax = axs[0]
    ax.plot(ht[:, 0], ht[:, 1], color=BLUE, label="head to tail")
    ax.plot(sbs[:, 0], sbs[:, 1], color=ORANGE, label="side by side")
    ax.axvspan(1000, 2000, color=BLUE, alpha=0.07, lw=0)
    ax.text(1500, 11.3, "head-to-tail data\nnot used here", ha="center", fontsize=8, color=BLUE)
    i = np.argmin(ht[:, 1])
    ax.annotate(f"closest: {ht[i, 1]:.1f} Å\n(atoms 2.3 Å apart)", xy=(ht[i, 0], ht[i, 1]),
                xytext=(ht[i, 0] - 650, 3.6), fontsize=8, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1))
    ax.set_xlabel("time (fs)")
    ax.set_ylabel("distance between centres (Å)")
    ax.set_ylim(2.5, 12.5)
    ax.set_title("Separation of the two waters", fontsize=10)
    ax.legend(fontsize=9, loc="lower left")
    ax = axs[1]
    ax.plot(sbs[:, 0], sbs[:, 2], color=ORANGE, label="water 1")
    ax.plot(sbs[:, 0], sbs[:, 3], color=ORANGE, ls="--", label="water 2")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("time (fs)")
    ax.set_ylabel("turn away from cavity axis (deg)")
    ax.set_title("Side by side: the waters swing as mirror images", fontsize=10)
    ax.set_ylim(-27, 27)
    ax.legend(fontsize=9, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "motion.png", dpi=200)
    plt.close(fig)


def ratio_rows(pk):
    def r(a, b, sa, sb):
        x = pk[a]["rabi"] / pk[b]["rabi"]
        # each peak uncertain by +-s cm^-1 -> Rabi by +-sqrt(2) s
        e = x * np.hypot(np.sqrt(2) * sa / pk[a]["rabi"], np.sqrt(2) * sb / pk[b]["rabi"])
        return x, e
    return [
        ("q = 0 start\nhead to tail, 0.4 ps", *r("ht_q0", "1w_q0", 10, 10), GREY),
        ("relaxed\nhead to tail, 0.4 ps", *r("ht_rel", "1w_rel", 10, 10), BLUE),
        ("relaxed\nside by side, 0.4 ps", *r("sbs_rel", "1w_rel", 10, 10), ORANGE),
        ("relaxed\nhead to tail, first 1 ps", *r("ht_rel_long_1ps", "1w_rel_long_1ps", 7, 7), BLUE),
        ("relaxed\nside by side, 2 ps", *r("sbs_rel_long", "1w_rel_long", 5, 5), ORANGE),
    ]


def plot_ratio(rows):
    fig, ax = plt.subplots(figsize=(9, 3.8))
    for i, (lab, x, e, c) in enumerate(rows):
        ax.errorbar(i, x, yerr=e, fmt="o", color=c, ms=8, capsize=4, lw=1.5)
        ax.text(i + 0.12, x, f"{x:.2f}", fontsize=9, va="center", color=c)
    ax.axhline(np.sqrt(2), color=INK, ls="--", lw=1)
    ax.text(0.6, np.sqrt(2) + 0.012, "√2 = 1.41 (expected)", fontsize=9, ha="left")
    ax.axvline(0.5, color=INK2, lw=0.8, ls=":")
    ax.set_xticks(range(len(rows)), [r[0] for r in rows], fontsize=8)
    ax.set_ylabel("Rabi(2 waters) / Rabi(1 water)")
    ax.set_xlim(-0.5, len(rows) - 0.4)
    ax.set_ylim(1.1, 1.75)
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(OUT / "ratio.png", dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pk = {k: peaks(p) for k, p in RUNS.items()}
    pk["ht_rel_long_1ps"] = peaks(RUNS["ht_rel_long"], FIRST_PS + 1)
    pk["1w_rel_long_1ps"] = peaks(RUNS["1w_rel_long"], FIRST_PS + 1)
    with open(OUT / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["key", "run", "steps_used", "resolution_cm-1", "LP_or_bend_cm-1", "UP_cm-1", "Rabi_cm-1"])
        for k, p in pk.items():
            run = RUNS[k.replace("_1ps", "")].relative_to(REPO / "runs")
            w.writerow([k, run, p["steps"], f"{p['res']:.1f}", f"{p['lp']:.0f}",
                        f"{p['up']:.0f}" if p["rabi"] else "", f"{p['rabi']:.0f}" if p["rabi"] else ""])
            print(f"{k:18s} {p['steps']:5d}  LP {p['lp']:6.0f}  UP {p['up']:6.0f}  Rabi {p['rabi']}")
    rows = ratio_rows(pk)
    for lab, x, e, _ in rows:
        print(f"{lab.replace(chr(10), ' '):40s} {x:.2f} +- {e:.2f}")
    plot_arrangements()
    plot_photon_start(pk)
    plot_spectra(pk)
    plot_motion(motion())
    plot_ratio(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
