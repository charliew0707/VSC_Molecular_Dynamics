"""
Plot a CO2 lambda-sweep campaign: IR spectra, polariton branches, Rabi splitting.

Reads every run folder in a campaign (made by scripts/campaigns/co2_lambda_sweep.sh),
groups them into series by driver + chi tag (folder names *_chi / *_nochi / *_bare),
and compares against the group's benchmark in polaritonic_deep_md
(tests/e2e/compare_with_bonini2024/data: NEP chi included/neglected + octopus QEDFT).

Peaks are found the same way as the group benchmark (plot.polariton_peaks): FFT of
the mean-subtracted dipole_x, two most intense local maxima in 1200-3300 cm^-1.

Usage:
    python scripts/plot_co2_sweep.py runs/co2/2026-10-02_lambda_sweep
Writes to results/co2/<campaign>/:
    ir_spectra.png        stacked IR spectra vs lambda, one panel per series
    polariton_branches.png  LP / UP vs lambda, with group references
    rabi_splitting.png    UP - LP vs lambda, with group references
    peaks.csv             every run's LP, UP, Rabi splitting
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
SHARED = REPO.parent / "polaritonic_deep_md"
sys.path.insert(0, str(SHARED))
import scripts.infrared as infrared  # noqa: E402  (shared repo, read-only)

BENCH = SHARED / "tests/e2e/compare_with_bonini2024/data"
BAND = (1200.0, 3300.0)        # same band as the group benchmark
OMEGA_C = 2400.0

# Validated categorical palette (dataviz default, light mode), fixed order per series
SERIES = {
    "nep_chi":     dict(label="NEP, χ included",    color="#2a78d6", marker="o"),
    "nep_nochi":   dict(label="NEP, χ neglected",   color="#eb6834", marker="s"),
    "pyscf_nochi": dict(label="PySCF, χ neglected", color="#1baf7a", marker="^"),
}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "font.size": 10, "legend.frameon": False,
})


# ── data ──────────────────────────────────────────────────────────────────────

def spectrum(dipole_x, dt):
    w, ft = infrared.fourier_transform(dipole_x - dipole_x.mean(), dt)
    return w * infrared.P_cm1, np.abs(ft)


def polariton_peaks(freq, amp):
    m = (freq > BAND[0]) & (freq < BAND[1])
    f, a = freq[m], amp[m]
    loc = [k for k in range(1, len(a) - 1) if a[k] > a[k - 1] and a[k] > a[k + 1]]
    if not loc:
        loc = [int(np.argmax(a))]
    top = sorted(sorted(loc, key=lambda k: -a[k])[:2], key=lambda k: f[k])
    return f[top[0]], f[top[-1]]


def load_campaign(campaign):
    """{series: [run dicts sorted by lambda]}; bare runs go in as lambda = 0."""
    runs = {}
    def dipole_file(d):
        return d / "dipole.dat" if (d / "dipole.dat").exists() else d / "dipole.dat.gz"

    for d in sorted(p for p in campaign.iterdir() if dipole_file(p).exists()):
        cfg = json.loads((d / "in.json").read_text())
        n_expected = cfg["steps"] + 1
        dip = np.loadtxt(dipole_file(d))
        if len(dip) < n_expected:
            print(f"skipping {d.name}: incomplete ({len(dip)}/{n_expected} steps)")
            continue
        driver = cfg["driver"]
        lam = cfg["lambda_photon"][0] if cfg.get("photons") else 0.0
        freq, amp = spectrum(dip[:, 2], dip[1, 1] - dip[0, 1])
        if lam == 0:      # no cavity: one peak (the bare stretch), as in the benchmark
            m = (freq > BAND[0]) & (freq < BAND[1])
            lp = up = freq[m][np.argmax(amp[m])]
        else:
            lp, up = polariton_peaks(freq, amp)
        run = dict(name=d.name, driver=driver, lam=lam, steps=cfg["steps"],
                   freq=freq, amp=amp, lp=lp, up=up, res=freq[1] - freq[0])
        tags = ["chi", "nochi"] if d.name.endswith("_bare") else [d.name.rsplit("_", 1)[1]]
        for tag in tags:          # the bare run is the lambda = 0 point of every series
            key = f"{driver}_{tag}"
            if key in SERIES and (tag != "chi" or driver != "pyscf"):
                runs.setdefault(key, []).append(run)
    return {k: sorted(v, key=lambda r: r["lam"]) for k, v in runs.items()}


def load_reference():
    ref = {}
    for key, fname in [("nep_chi", "nep_chi_included.txt"), ("nep_nochi", "nep_chi_neglected.txt")]:
        a = np.loadtxt(BENCH / fname)
        ref[key] = dict(lam=a[0], lp=a[1], up=a[2])
    o = np.loadtxt(BENCH / "octopus_CO2_CBOA_FD.txt")   # rows: lambda, bend, LP, UP
    ref["qedft"] = dict(lam=o[0], lp=o[2], up=o[3])
    return ref


# ── figures ───────────────────────────────────────────────────────────────────

def plot_spectra(series, out):
    keys = [k for k in SERIES if k in series]
    fig, axes = plt.subplots(1, len(keys), figsize=(4.2 * len(keys), 6.5), sharey=True, squeeze=False)
    for ax, key in zip(axes[0], keys):
        s = SERIES[key]
        for i, r in enumerate(series[key]):
            m = (r["freq"] > 1000) & (r["freq"] < 3400)
            f, a = r["freq"][m], r["amp"][m]
            band = (f > BAND[0]) & (f < BAND[1])
            y = a / a[band].max() * 0.85 + i
            ax.plot(f, y, color=s["color"] if r["lam"] > 0 else INK2, lw=1.5)
            ax.text(1020, i + 0.12, "no cavity" if r["lam"] == 0 else f"λ = {r['lam']:g}",
                    fontsize=9, color=INK)
            if r["lam"] > 0:
                for pk in (r["lp"], r["up"]):
                    ax.plot([pk, pk], [i, i + 0.92], color=INK2, lw=0.8, ls=":")
        ax.axvline(OMEGA_C, color=INK2, lw=1, ls="--")
        ax.text(OMEGA_C + 20, 1.0, r"$\omega_c$", color=INK2, fontsize=10,
                transform=ax.get_xaxis_transform(), va="top")
        ax.set_title(s["label"], fontsize=11, loc="left")
        ax.set_xlabel("Frequency (cm⁻¹)")
        ax.set_xlim(1000, 3400)
        ax.set_yticks([])
        ax.grid(axis="y", visible=False)
    axes[0][0].set_ylabel("IR intensity (normalized, offset by λ)")
    fig.suptitle("CO₂ IR spectra in a cavity resonant with the asymmetric stretch "
                 rf"($\omega_c$ = {OMEGA_C:g} cm⁻¹); dotted lines = LP / UP", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out / "ir_spectra.png", dpi=160)
    plt.close(fig)


def _reference_lines(ax, ref, value):
    ax.plot(ref["qedft"]["lam"], value(ref["qedft"]), color=INK, lw=1.2, ls="-",
            label="QEDFT reference (octopus)", zorder=1)
    for key in ("nep_chi", "nep_nochi"):
        ax.plot(ref[key]["lam"], value(ref[key]), color=SERIES[key]["color"], lw=1.2, ls="--",
                alpha=0.8, label=f"group benchmark: {SERIES[key]['label']}", zorder=1)


def _points(ax, series, value, with_label=True):
    for key, s in SERIES.items():
        if key not in series:
            continue
        lam = [r["lam"] for r in series[key]]
        ax.plot(lam, [value(r) for r in series[key]], ls="none", marker=s["marker"], ms=8,
                color=s["color"], mec=SURFACE, mew=2, label=f"this campaign: {s['label']}" if with_label else None,
                zorder=3)


def plot_branches(series, ref, out):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    _reference_lines(ax, ref, lambda d: d["lp"])
    _reference_lines(ax, {k: dict(lam=v["lam"], lp=v["up"]) for k, v in ref.items()}, lambda d: d["lp"])
    _points(ax, series, lambda r: r["lp"])
    _points(ax, series, lambda r: r["up"], with_label=False)
    ax.axhline(OMEGA_C, color=INK2, lw=0.8, ls=":")
    ax.text(0.316, OMEGA_C, r"$\omega_c$", color=INK2, fontsize=10, va="center")
    ax.text(0.305, ref["qedft"]["up"][-1], "UP", color=INK2, fontsize=10, va="center")
    ax.text(0.305, ref["qedft"]["lp"][-1], "LP", color=INK2, fontsize=10, va="center")
    handles, labels = ax.get_legend_handles_labels()
    keep = {l: h for h, l in zip(handles, labels)}       # drop duplicate reference labels
    ax.legend(keep.values(), keep.keys(), fontsize=8.5, loc="lower left")
    ax.set_xlabel("Coupling strength λ (a.u.)")
    ax.set_ylabel("Polariton frequency (cm⁻¹)")
    ax.set_title("CO₂ polariton branches vs coupling", loc="left", fontsize=11)
    ax.set_xlim(-0.005, 0.33)
    fig.tight_layout()
    fig.savefig(out / "polariton_branches.png", dpi=160)
    plt.close(fig)


def plot_rabi(series, ref, out):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    _reference_lines(ax, ref, lambda d: d["up"] - d["lp"])
    _points(ax, series, lambda r: r["up"] - r["lp"])
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_xlabel("Coupling strength λ (a.u.)")
    ax.set_ylabel(r"Rabi splitting $\Omega_R$ = UP − LP (cm⁻¹)")
    ax.set_title("CO₂ Rabi splitting vs coupling", loc="left", fontsize=11)
    ax.set_xlim(-0.005, 0.32)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out / "rabi_splitting.png", dpi=160)
    plt.close(fig)


def write_table(series, out):
    with open(out / "peaks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "run", "lambda", "steps", "resolution_cm-1", "LP_cm-1", "UP_cm-1", "Rabi_cm-1"])
        for key, runs in series.items():
            for r in runs:
                w.writerow([key, r["name"], r["lam"], r["steps"], f"{r['res']:.1f}",
                            f"{r['lp']:.0f}", f"{r['up']:.0f}", f"{r['up'] - r['lp']:.0f}"])


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    campaign = Path(sys.argv[1]).resolve()
    out = REPO / "results" / "co2" / campaign.name
    out.mkdir(parents=True, exist_ok=True)
    series, ref = load_campaign(campaign), load_reference()
    plot_spectra(series, out)
    plot_branches(series, ref, out)
    plot_rabi(series, ref, out)
    write_table(series, out)
    for key, runs in series.items():
        print(f"\n{SERIES[key]['label']}  ({runs[-1]['steps']} steps, {runs[-1]['res']:.1f} cm⁻¹ resolution)")
        for r in runs:
            print(f"  λ = {r['lam']:<5g}  LP {r['lp']:6.0f}  UP {r['up']:6.0f}  Ω_R {r['up'] - r['lp']:5.0f} cm⁻¹")
    print(f"\nwrote {out.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
