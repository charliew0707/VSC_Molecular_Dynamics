# Running VSC simulations: a guide

How to run cavity MD with the Flick group's `polaritonic_deep_md` code from this
repo: where everything lives, where results are saved, and what to watch out for.

---

## 1. The big picture

```
~/Flick Group Code/
├── polaritonic_deep_md/          ← SHARED group code. Read-only: only ever `git pull`.
└── VSC_Molecular_Dynamics/       ← YOUR repo: inputs, runs, data, plots, reports.
    ├── configs/                  ← presets that run.py combines into in.json
    │   ├── defaults.json         ←   steps, timestep, nphoton
    │   ├── molecules/co2.json    ←   geometry, cavity frequency (cm⁻¹), polarization
    │   ├── molecules/h2o.json
    │   └── drivers/{pyscf,nep,ff}.json  ← backend settings (xc, basis, polar...)
    ├── geometries/               ← .xyz starting structures
    ├── runs/<molecule>/<run>/    ← ★ every run's inputs + raw output data land here
    ├── scripts/run.py            ← the launcher
    ├── scripts/*.py              ← analysis / plotting
    ├── results/<molecule>/       ← finished figures you want to keep
    ├── reports/                  ← write-ups (.tex, .pdf)
    ├── large/                    ← big files git ignores (NEP models, huge trajectories)
    └── environment.yml           ← the `vsc` conda environment
```

**Rule of thumb:** the group code does the physics. Everything you create goes
in `VSC_Molecular_Dynamics`.

---

## 2. One-time setup (already done on this laptop)

```bash
cd ~/"Flick Group Code/VSC_Molecular_Dynamics"
conda env create -f environment.yml
conda activate vsc
pip install -e ../polaritonic_deep_md
```

`pip install -e` *links* the shared code into the environment instead of copying
it, so a `git pull` there takes effect immediately with no reinstall.

Don't use the base Anaconda environment. Its scipy is broken and `cboamd`
crashes on import.

---

## 3. Doing a run

Start every new terminal with:

```bash
conda activate vsc
cd ~/"Flick Group Code/VSC_Molecular_Dynamics"
```

Then:

```bash
python scripts/run.py <molecule> <driver> [options]
```

### Common recipes

| Goal | Command |
|---|---|
| Bare molecule (no cavity), reference spectrum | `python scripts/run.py h2o pyscf` |
| Molecule in a cavity on resonance | `python scripts/run.py h2o pyscf --cavity --lam 0.1` |
| Coupling-strength sweep (one run per λ) | `python scripts/run.py h2o pyscf --cavity --lam 0.02 0.05 0.1 0.2` |
| Detuned cavity | `python scripts/run.py co2 ff --cavity --lam 0.1 --omega 2300` |
| Longer run, better spectral resolution | `python scripts/run.py h2o pyscf --cavity --lam 0.1 --steps 8000` |
| Different DFT basis/functional | `python scripts/run.py h2o pyscf --set basis=augccpvdz --set xc=pbe` |
| Label a run | `... --name first_try` |
| Check the input without running | `... --dry-run` |
| Skip the IR spectrum step | `... --no-ir` |

### All options

| Option | Meaning | Default |
|---|---|---|
| `molecule` | a file in `configs/molecules/` (`co2`, `h2o`) | required |
| `driver` | a file in `configs/drivers/` (`pyscf`, `nep`, `ff`) | required |
| `--cavity` | turn photons on | off |
| `--lam L [L ...]` | coupling strength(s) λ in a.u.; several values → several runs | `0.05` |
| `--omega W` | cavity frequency in **cm⁻¹** | molecule preset (CO₂ 2430, H₂O 1584) |
| `--steps N` | MD steps | `2000` |
| `--dt T` | timestep in fs | `0.5` |
| `--set KEY=VALUE` | override any `in.json` key (repeatable) | none |
| `--name TAG` | append a tag to the folder name | none |
| `--dry-run` | print `in.json` and folder name only | off |
| `--no-ir` | don't run `infrared` afterwards | off |

`--set` accepts JSON values: `--set steps=500`, `--set polar=true`,
`--set 'lambda_vector=[[0,1,0]]'`.

### Drivers

| Driver | What it is | Speed | Works for |
|---|---|---|---|
| `ff` | built-in classical CO₂ force field | seconds | **CO₂ only** |
| `pyscf` | DFT on the fly (LDA / cc-pVDZ by default) | minutes to hours | any molecule |
| `nep` | trained neural-network potential | fast (200 cavity steps ≈ 1 s) | only where a model exists: CO₂ (see §7) |

**Rough cost:** 10 PySCF H₂O steps take about 3 s on this laptop, so the default
2000 steps takes about 10 minutes. CO₂ `ff` runs 200 steps in well under a second.

---

## 4. Where the data goes

Each run creates **its own new folder**:

```
runs/<molecule>/<date>_<driver>[_cav_w<ω cm⁻¹>_lam<λ>][_<name>]/
```

For example:

```
runs/h2o/2026-09-29_pyscf/                         bare H₂O, PySCF
runs/h2o/2026-09-29_pyscf_cav_w1584_lam0.1/        H₂O in cavity, ω=1584 cm⁻¹, λ=0.1
runs/co2/2026-09-29_ff_cav_w2430_lam0.05_sweep1/   tagged with --name sweep1
```

The launcher **refuses to overwrite** an existing folder. To repeat an identical
run on the same day, give it a `--name`.

### What's in a run folder

| File | Contents | Units |
|---|---|---|
| `in.json` | exact input given to `cboamd` | |
| `run_info.json` | command used, start time, machine, Python, **shared-code commit** | |
| `run.log` | everything printed during the run (check here if it fails) | |
| `settings.json` | the full settings `cboamd` actually used, including defaults | |
| `md.log` | per-step time, total/potential/kinetic energy, temperature | ps, eV, K |
| `dipole.dat` | 8 columns (header lists only 5): 0 step, 1 time, **2–4 dipole x/y/z used by `infrared`** (`dipolepol` in the code), 5–7 a second dipole x/y/z (`dipole`; `compare_bend_hanning.py` uses col 7) | a.u. |
| `energy.dat` | step, time, energy | eV |
| `photon.dat` | step, time, then qa, pa, Ea for each cavity mode | a.u. |
| `polarizability.dat` | step, time, 3×3 polarizability tensor | a.u. |
| `force.dat` | step, time, force on every atom (x,y,z) | eV/Å |
| `force_bare.dat` | same but without the cavity field (cavity runs only) | eV/Å |
| `position.dat` | step, time, coordinates (+ velocities) | Å |
| `md.traj` | ASE trajectory: open with `ase gui md.traj` or `ase.io.read` | |
| `velocities.xyz` | final velocities | |
| `spectrum_ase.dat` | **IR spectrum**: col 1 = frequency (cm⁻¹), col 2 = intensity | cm⁻¹, arb. |

Load any `.dat` with `np.loadtxt("dipole.dat")`. Lines starting with `#` are
skipped automatically.

---

## 5. Analysis and figures

- Write analysis scripts in `scripts/` and point them at run folders. For example,
  `scripts/compare_bend_hanning.py` compares a bare and a cavity spectrum:
  ```bash
  cd runs/h2o && python ../../scripts/compare_bend_hanning.py
  ```
- Save figures you want to keep in `results/<molecule>/`, and reference them from
  `reports/vsc_report.tex`. It already looks in `results/co2/` and `results/h2o/`.
- For notebooks, run `jupyter lab` inside the `vsc` environment. In VS Code,
  select the **vsc** interpreter.

### Physics tips

- **Rabi splitting** Ω_R = UP − LP. Compare a bare run with a cavity run at the
  same settings and find the two peaks either side of the bare peak.
- **Spectral resolution** is about 33,356 / (steps × dt [fs]) cm⁻¹. The default
  2000 × 0.5 fs gives about 33 cm⁻¹, which is fine for H₂O's ~125 cm⁻¹ splitting.
  For small splittings (small λ), use `--steps 8000` or more.
- **`infrared -m ase` only uses the x-component of the dipole** (`dipole.dat`
  column 3). So orient molecules so the mode you care about has its dipole
  along **x**, and polarize the cavity along x (`lambda_vector: [[1,0,0]]`). Both
  presets do this: CO₂ lies along x, and H₂O (`geometries/h2o-single-x.xyz`) has
  its symmetry axis along x, which is the bend's dipole direction.
  - The older H₂O runs (`runs/h2o/pyscf`, `runs/h2o/pyscf_cavity`) used the
    z-oriented `h2o-single.xyz`. Their `spectrum_ase.dat` files have been
    recomputed from the z column (see the header line). The Hanning analysis
    (`compare_bend_hanning.py`, 125 cm⁻¹ splitting) always read z directly and
    was unaffected.
  - Quick check for any run:
    `python -c "import numpy as np; print(np.loadtxt('dipole.dat')[:,2:5].std(0))"`.
    The largest value should be x.
- Converting cavity frequency: ω[a.u.] = ω[cm⁻¹] / 219474.63. `--omega` does this
  for you.

---

## 6. Adding things

- **A new molecule:** put the `.xyz` in `geometries/`, then create
  `configs/molecules/<name>.json`:
  ```json
  {
      "xyz_file": "geometries/<name>.xyz",
      "omega_cm": 1650,
      "lambda_vector": [[0, 0, 1]]
  }
  ```
  After that, `python scripts/run.py <name> pyscf ...` works.
- **A new driver preset,** such as a PBE variant: copy `configs/drivers/pyscf.json`
  to `configs/drivers/pyscf_pbe.json` and edit it. It then shows up as the
  driver `pyscf_pbe`.
- **Molecule-specific driver settings** go inside the molecule file under the
  driver's name. See the `"nep"` block in `co2.json`.
- **Changing defaults for everything:** edit `configs/defaults.json`.

---

## 7. NEP models

The `nep` driver for CO₂ uses the NEP models (energy, dipole, polarizability)
that Johannes committed to the shared repo on 2026-06-15:

```
../polaritonic_deep_md/tests/e2e/no-polar/mlip/models/nep-{energy,dipole,polar}.txt
```

`configs/molecules/co2.json` points straight at them, so a `git pull` there
picks up any updated models. They reproduce the older `runs/co2/nep` run exactly.
They're the models the group's own tests and Bonini-comparison benchmarks use.
Yli's production models (`/mnt/home/yli11/...` on the cluster) may differ. To
use other models, drop them in `large/models/` and override with
`--set nep_pot=... --set nep_dip=... --set nep_pol=...`.

NEP runs print an ASE `FutureWarning` about `ignore_bad_restart_file`. It's
harmless.

---

## 8. Keeping the shared code up to date (and safe)

Update:

```bash
cd ~/"Flick Group Code/polaritonic_deep_md" && git pull --ff-only
```

- The shared repo has a **local lock**: commits and pushes are blocked on this
  machine, but pulling works.
  - Override once: `git commit --no-verify ...`
  - Remove: `rm .git/hooks/pre-commit .git/hooks/pre-push` in that repo
- **Never edit or run inside the shared repo.** `cboamd` writes its outputs into
  whatever folder it's run from. `run.py` always runs in a folder under `runs/`.
- Every run records the shared-code commit in `run_info.json`. If results change
  after a `git pull`, compare commits:
  `git -C ../polaritonic_deep_md log --oneline <old>..<new>`.

---

## 9. Saving your work (git)

```bash
cd ~/"Flick Group Code/VSC_Molecular_Dynamics"
git add runs/h2o/2026-09-29_pyscf_cav_w1584_lam0.1 results/h2o/new_plot.png
git commit -m "H2O cavity run lambda=0.1"
git push
```

- **Commit:** run folders with small outputs, figures, scripts, configs and reports.
- **Don't commit** anything over about 50 MB (GitHub blocks files over 100 MB).
  Move big trajectories to `large/` and keep only the `.dat` files and
  `spectrum_ase.dat` in the run folder.
- Check sizes before committing with `du -sh runs/*/*`.
- Delete failed or test runs before committing: `rm -rf runs/co2/<bad-run>`.

---

## 10. Troubleshooting

| Symptom | Fix |
|---|---|
| `ImportError ... scipy` / `_spropack` | You're in base Anaconda. Run `conda activate vsc`. |
| `'cboamd' not found next to ...python` | Wrong environment, or the group code isn't installed in it: `pip install -e ../polaritonic_deep_md` |
| `... already exists; use --name` | Same run already done today. Add `--name v2`. |
| `warning: nep_pot not found` | Model path is wrong; check that the shared repo is at `../polaritonic_deep_md` (§7). |
| Run crashed partway | Read `run.log` in the run folder, then delete the folder and re-run. |
| `polaritonic_deep_md is locked` | You tried to commit in the shared repo. Commit in your own repo instead. |
| Claude can't see a folder | macOS blocks Desktop, Documents and Downloads. Keep projects under `~/Flick Group Code/`. |
| Results differ after updating shared code | Compare the `polaritonic_deep_md.commit` in each run's `run_info.json`. |
