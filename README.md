# VSC Molecular Dynamics

Cavity Born–Oppenheimer MD simulations of molecules inside optical microcavities
(vibrational strong coupling, VSC). This is Charlie Wade's personal repo for runs,
data, analysis, and write-ups. The simulations use the Flick group's shared code,
[`polaritonic_deep_md`](https://github.com/flickgroup/polaritonic_deep_md).

The shared code is **not** copied or edited here. It lives next to this repo:

```
Flick Group Code/
├── polaritonic_deep_md/      # shared group repo (read-only for me, just `git pull`)
└── VSC_Molecular_Dynamics/   # this repo
```

## Results

| Molecule | Mode | Rabi splitting (λ=0.1) |
|----------|------|------------------------|
| CO₂ | Asymmetric C–O stretch (~2349 cm⁻¹) | ~420 cm⁻¹ |
| H₂O | Bend (~1580 cm⁻¹) | 125 cm⁻¹ |

## Layout

| Folder        | What goes in it |
|---------------|-----------------|
| `geometries/` | Input structures (`.xyz`) |
| `runs/`       | One folder per run: `runs/<molecule>/<run-name>/` with its `in.json` and outputs |
| `configs/`    | Run presets (defaults, molecules, drivers) used by `scripts/run.py` |
| `scripts/`    | `run.py` launcher, analysis and plotting scripts |
| `results/`    | Figures produced from runs |
| `reports/`    | Write-ups (`vsc_report.tex`, PDFs) |
| `docs/`       | Guide and notes |
| `large/`      | Big trajectories/checkpoints, which git ignores |

### Runs so far

- `runs/co2/{ff,nep,pyscf}`: CO₂ force-backend comparison (FF vs PySCF vs NEP dipoles)
- `runs/co2/2026-06-09_pyscf_test`, `runs/co2/co2_test_001`: early PySCF cavity test runs
- `runs/h2o/pyscf`, `runs/h2o/pyscf_cavity`: H₂O bare vs cavity spectra (bend Rabi splitting)

**Full guide: [docs/GUIDE.md](docs/GUIDE.md)**: running, where data goes, analysis tips, troubleshooting.

## Setup (once)

```bash
conda env create -f environment.yml
conda activate vsc
pip install -e ../polaritonic_deep_md
```

After that, run `conda activate vsc` in each new terminal.

## Running

`scripts/run.py` builds `in.json` from the presets in `configs/`, makes a new
dated folder in `runs/<molecule>/`, and runs `cboamd` and then `infrared` there:

```bash
python scripts/run.py h2o pyscf                            # bare molecule
python scripts/run.py h2o pyscf --cavity --lam 0.1         # cavity on resonance
python scripts/run.py co2 nep --cavity --lam 0.05 0.1 0.2  # coupling sweep
python scripts/run.py co2 ff --steps 20 --dry-run          # just print in.json
python scripts/run.py co2 pyscf --set basis=augccpvdz --name bigbasis
python scripts/run.py --help
```

Each run folder gets `in.json`, `run_info.json` (command + shared-repo commit),
`run.log`, and all cboamd outputs (`dipole.dat`, `spectrum_ase.dat`, ...).

- **Presets:** `configs/defaults.json`, `configs/molecules/<mol>.json` (geometry,
  cavity frequency in cm⁻¹, polarization), `configs/drivers/<driver>.json`. Add a
  JSON file there to add a molecule or driver.
- **NEP models** are not in git. Put the CO₂ models in `large/models/co2_nep/`
  (see `configs/molecules/co2.json` for the paths).

## Analysis

```bash
cd runs/h2o
python ../../scripts/compare_bend_hanning.py
```

Produces `h2o_bend_hanning.png` and `h2o_full_hanning.png`.

## Updating the shared code

```bash
cd ~/"Flick Group Code/polaritonic_deep_md" && git pull --ff-only
```

That repo has a local lock (git hooks) that blocks commits and pushes from this
machine. Pulls are unaffected. To unlock it, delete `.git/hooks/pre-commit` and
`.git/hooks/pre-push` in that repo.

## Reproducibility

Shared code was at commit `f0c0ed7` (2026-07-24) when this README was written.
Runs before that date used older versions. When starting a new run, note the
shared-repo commit (`git -C ../polaritonic_deep_md rev-parse --short HEAD`) in the
run folder.
