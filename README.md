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
| `scripts/`    | Analysis and plotting scripts |
| `results/`    | Figures produced from runs |
| `reports/`    | Write-ups (`vsc_report.tex`, PDFs) |
| `docs/`       | Notes |
| `large/`      | Big trajectories/checkpoints, which git ignores |

### Runs so far

- `runs/co2/{ff,nep,pyscf}`: CO₂ force-backend comparison (FF vs PySCF vs NEP dipoles)
- `runs/co2/2026-06-09_pyscf_test`, `runs/co2/co2_test_001`: early PySCF cavity test runs
- `runs/h2o/pyscf`, `runs/h2o/pyscf_cavity`: H₂O bare vs cavity spectra (bend Rabi splitting)

## Quick start

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
