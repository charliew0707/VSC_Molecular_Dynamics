# CLAUDE.md

Charlie Wade's personal research repo for vibrational strong coupling (VSC): cavity
Born–Oppenheimer MD of molecules in an optical cavity. It holds **runs, data, plots and
reports**. The physics code is the Flick group's shared repo, which lives next to this one
and is only *used* from here.

```
~/Flick Group Code/
├── polaritonic_deep_md/      # shared group code (GitHub flickgroup), READ-ONLY
└── VSC_Molecular_Dynamics/   # this repo (GitHub charliew0707)
```

`README.md` has the results and `docs/GUIDE.md` is the full user guide. Keep both
updated when you add a campaign or change the workflow.

## Hard rules

- **Never edit, commit to, or run inside `../polaritonic_deep_md`.** The only allowed
  command there is `git pull --ff-only`. It has local pre-commit/pre-push hooks that
  block commits on purpose; don't remove or bypass them. Read its code freely.
- **Stage explicit paths, never `git add -A` / `git add .`.** Charlie runs his own tests
  in this repo, and a blanket add once pushed a half-finished run of his. Check
  `git status` and ask about run folders you didn't create.
- **Ask before deleting run data** you didn't create in this session.
- Use the **`vsc` conda env** (`/opt/anaconda3/envs/vsc/bin/python`). The group code is
  `pip install -e` linked there. Don't install into or run from base Anaconda.
- **Background jobs started by Claude are killed after ~30 min.** Size every batch to
  finish in under ~25 min (see costs below). Longer runs: give Charlie the command to run
  in his own terminal.

## Running

```bash
python scripts/run.py <molecule> <driver> [--cavity] [--lam L ...] [--omega CM] \
    [--steps N] [--set key=value] [--campaign NAME] [--name TAG] [--dry-run]
```

- Builds `in.json` from `configs/defaults.json` + `configs/molecules/<mol>.json` +
  `configs/drivers/<driver>.json` + CLI options. Makes a fresh folder
  `runs/<mol>/[<campaign>/]<run>/`, writes `run_info.json` (command + shared-repo commit),
  runs `cboamd` then `infrared -m ase`, and writes `dipole.dat.gz` / `energy.dat.gz`.
  It refuses to overwrite an existing folder.
- A sweep is a **campaign**: `scripts/campaigns/<name>.sh` (re-creates the runs) +
  `scripts/plot_<name>.py` (writes figures and `peaks.csv` to `results/<mol>/<campaign>/`).
  Campaign folders are date-prefixed, e.g. `2026-10-02_lambda_sweep`.
- Run folder suffixes the plot scripts rely on: `_bare` (no cavity, λ = 0), `_chi`
  (polarizability on), `_nochi` (off).
- Parallel PySCF: one process per run with `OMP_NUM_THREADS=1`.

### Approximate costs (8-core laptop)

| System | Driver | Per step, alone | Notes |
|---|---|---|---|
| CO₂ | NEP | ~5 ms | 10,000 steps ≈ 1 min; the default for CO₂ |
| CO₂ | PySCF, χ off | ~1.1 s | 5 runs × 800 steps in parallel ≈ 23 min (the limit) |
| H₂O | PySCF, χ off | ~0.17 s | 7 runs × 800 steps in parallel ≈ 7 min |
| H₂O | PySCF, χ on | ~0.8 s | batches of ≤ 3 × 800 steps ≈ 15 min |

Frequency resolution = 33356 / (steps × dt[fs]) cm⁻¹: 800 steps → 83, 10,000 → 6.7.

## Git / data policy

`.gitignore` keeps only `in.json`, `settings.json`, `run_info.json`, `dipole.dat.gz` and
`energy.dat.gz` from each run (enough to regenerate every spectrum and plot). Full
`.dat`, `md.traj` and logs stay on disk only. Plot scripts read `dipole.dat` if present,
else `dipole.dat.gz`. Runs committed before 2026-10-02 still have older tracked files.

## Physics / code gotchas (learned the hard way)

- **`infrared -m ase` uses only dipole x** (`dipole.dat` col 2). Orient the mode's
  dipole and the cavity polarization along x. CO₂ lies along x; H₂O uses
  `geometries/h2o-single-x.xyz` (C₂ axis along x). The old `h2o-single.xyz` is
  z-oriented, so its spectra came from noise.
- `dipole.dat` has 8 columns (the header lists 5): 2–4 = dipole used by `infrared`,
  5–7 = a second dipole.
- **`md.log` Etot excludes the cavity.** Conserved H = Etot + ½p² + ½(1+λ²χ)·ea²
  (p, ea from `photon.dat` cols 3–4). Use `scripts/check_energy.py`; ~1% fluctuation and
  < 0.5% drift is good. dt = 0.5 fs is fine.
- **χ (polarizability) on** screens the cavity to ω_c/√(1+λ²χ), red-detuning it. One
  polariton dominates and the other can be too weak to see in short runs. The group's
  "resonant" variant scales ω_c by √(1+λ²χ).
- **Dipole self-energy shifts the resonance:** CO₂ (NEP) bare stretch is 2438 but the
  effective one is 2464, so true resonance is ≈ 2500 cm⁻¹. CO₂ detuning matches the exact
  linear (Hopfield + self-energy) model to 2.6 cm⁻¹. H₂O sits 40–60 cm⁻¹ below it (open
  question).
- **Molecules rotate** when the polarization is tilted from their dipole (H₂O strongly,
  CO₂ beyond ~45°), so tilted-cavity runs are not clean cos θ tests.
- Short runs: unwindowed FFT side lobes are ~20% of the main peak. H₂O peak finding uses
  zero-padding plus a ≥ 50% height rule (`plot_h2o_sweep.py`). CO₂ uses the group
  benchmark's method (`plot_co2_sweep.py`).
- CO₂ NEP models = `../polaritonic_deep_md/tests/e2e/no-polar/mlip/models/`. NEP runs
  print a harmless ASE `FutureWarning`.

## Validation anchors

- CO₂ NEP λ sweep reproduces the group's `tests/e2e/compare_with_bonini2024` branches
  exactly (e.g. λ = 0.1, χ on: LP 2105 / UP 2525, Ω_R = 420 cm⁻¹).
- H₂O PySCF (LDA/cc-pVDZ, χ off) λ = 0.1: Ω_R = 125 cm⁻¹; Ω_R ≈ 1265·λ cm⁻¹.

## Plots and reports

- Matplotlib style: fixed categorical colours per series (NEP χ on `#2a78d6` ●,
  NEP χ off `#eb6834` ■, PySCF χ off `#1baf7a` ▲, PySCF χ on `#eda100` ●), shared
  `INK`/`SURFACE` constants in `plot_co2_sweep.py`. Look at each PNG before committing.
- Reports: LaTeX in `reports/`, built with `tectonic <file>.tex` (installed in `vsc`).
  Figures come from `../results/` via `\graphicspath`. Commit `.tex` + `.pdf`.
- Commit messages end with the `Co-Authored-By` line the harness provides.
