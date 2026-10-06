# VSC Molecular Dynamics

Charlie Wade's research repository for **vibrational strong coupling (VSC)**:
cavity Born–Oppenheimer molecular dynamics of molecules (CO₂, H₂O, …) coupled to an
optical cavity mode. It stores the inputs, raw data, analysis, figures and reports from
those simulations.

The simulations themselves use the Flick group's shared code,
[`polaritonic_deep_md`](https://github.com/flickgroup/polaritonic_deep_md). That code
is **not** copied or modified here. It sits next to this repo and is only used from it:

```
Flick Group Code/
├── polaritonic_deep_md/      # shared group code, read-only (just `git pull`)
└── VSC_Molecular_Dynamics/   # this repo
```

## How to read this repo

| Folder | What's in it |
|---|---|
| [`reports/`](reports/) | **Start here.** Write-ups of what was run and what it showed (LaTeX + PDF), named by date. |
| [`results/`](results/) | Figures and peak tables (`peaks.csv`), grouped by molecule, then by campaign. |
| [`runs/`](runs/) | Raw simulation data, grouped by molecule, then by campaign, one folder per run. |
| [`scripts/`](scripts/) | `run.py` (launches a run), `campaigns/` (scripts that re-create each set of runs), `plot_*.py` (turn runs into figures), `check_energy.py`. |
| [`configs/`](configs/) | Presets that `run.py` combines into each run's input: defaults, molecules, drivers. |
| [`geometries/`](geometries/) | Starting structures (`.xyz`). |
| [`docs/`](docs/) | [`GUIDE.md`](docs/GUIDE.md): how to set up, run, analyse and save everything. |

### Naming

- **Campaigns** are related sets of runs, in folders named `YYYY-MM-DD_<description>`
  (e.g. `runs/h2o/2026-10-02_lambda_sweep/`). Their figures go in the matching folder
  under `results/<molecule>/`. The script that re-creates a campaign is
  `scripts/campaigns/<molecule>_<description>.sh`, and its header names the
  `scripts/plot_*.py` that makes the figures.
- **Run folders** describe their settings:
  `<driver>[_cav_w<cavity cm⁻¹>_lam<coupling λ>][_<tag>]`. For example,
  `pyscf_cav_w1584_lam0.1_nochi` is PySCF, cavity at 1584 cm⁻¹, λ = 0.1, polarizability
  neglected. `_bare` means no cavity.
- Every run folder has `in.json` (exact input), `run_info.json` (command, date, and
  shared-code commit) and compressed `dipole.dat.gz` / `energy.dat.gz`. Bulky raw outputs
  stay on the machine that ran them.
- Undated folders directly under `runs/<molecule>/` are early runs from before the
  launcher existed.

## Using it

See **[docs/GUIDE.md](docs/GUIDE.md)** for environment setup, running simulations and
sweeps, where outputs go, making plots, energy checks, and troubleshooting.
