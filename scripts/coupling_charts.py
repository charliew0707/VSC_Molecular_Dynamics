"""
Polaritonic MD - Coupling Strength Comparison Plot
Generates a 3-panel figure comparing dipole_x, qa, and Ekin across runs.

Usage:
    python plot_polaritonic_md.py

Edit the RUNS dict below to point to your step summary CSV files.
"""

import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ── Configuration ─────────────────────────────────────────────────────────────

RUNS = {
    "λ = 0.01": "/Users/charleswade/Desktop/MD_coupling_runs_data/co2_0.01_step_summary.csv",
    "λ = 0.05": "/Users/charleswade/Desktop/MD_coupling_runs_data/md_step_summary_0.05.csv",
    "λ = 2.0":  "/Users/charleswade/Desktop/MD_coupling_runs_data/co2_2.0_step_summary.csv",
}

COLORS = {
    "λ = 0.01": "#60a5fa",
    "λ = 0.05": "#34d399",
    "λ = 2.0":  "#f472b6",
}

PANELS = [
    ("dipole_x_au", "Dipole Moment X (A.U.)"),
    ("qa",          "Photon Coordinate $q_a$"),
    ("Ekin_eV",     "Kinetic Energy $E_{kin}$ (eV)"),
]

# ── Load data ──────────────────────────────────────────────────────────────────

def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

data = {}
for label, path in RUNS.items():
    rows = load_csv(path)
    data[label] = {col: [float(r[col]) for r in rows] for col in rows[0].keys()}

# ── Plot ───────────────────────────────────────────────────────────────────────

plt.style.use("dark_background")

fig = plt.figure(figsize=(9, 10))
fig.patch.set_facecolor("#0f172a")

gs = gridspec.GridSpec(3, 1, hspace=0.45, left=0.12, right=0.95, top=0.93, bottom=0.07)

fig.suptitle(
    "Polaritonic MD — Coupling Strength Comparison\nCO₂  |  Steps 0–10  |  dt = 20.67 fs",
    fontsize=12, color="#e2e8f0", y=0.98, linespacing=1.6
)

for i, (col, ylabel) in enumerate(PANELS):
    ax = fig.add_subplot(gs[i])
    ax.set_facecolor("#1e293b")

    for label in RUNS:
        t = data[label]["time_fs"]
        y = data[label][col]
        ax.plot(t, y, color=COLORS[label], linewidth=2, marker="o",
                markersize=4, label=label, zorder=3)

    # zero line if data crosses zero
    yvals = np.concatenate([data[l][col] for l in RUNS])
    if yvals.min() < 0 < yvals.max():
        ax.axhline(0, color="white", linewidth=0.7, linestyle="--", alpha=0.3, zorder=2)

    ax.set_ylabel(ylabel, fontsize=10, color="#94a3b8", labelpad=8)
    ax.set_xlabel("Time (fs)", fontsize=9, color="#64748b")
    ax.tick_params(colors="#94a3b8", labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor("#334155")
    ax.grid(color="#334155", linewidth=0.5, zorder=1)

    if i == 0:
        ax.legend(
            loc="upper right", fontsize=9,
            framealpha=0.3, edgecolor="#334155",
            labelcolor="#e2e8f0"
        )

plt.savefig("polaritonic_md_comparison.png", dpi=150, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print("Saved: polaritonic_md_comparison.png")
plt.show()