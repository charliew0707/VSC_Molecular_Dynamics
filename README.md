# VSC_Molecular_Dynamics

Cavity Born–Oppenheimer MD simulations of molecules inside optical microcavities,
using the [polaritonic_deep_md](https://github.com/flickgroup) code (Flick group).

## Results

| Molecule | Mode | Rabi splitting (λ=0.1) |
|----------|------|------------------------|
| CO₂ | Asymmetric C–O stretch (~2349 cm⁻¹) | ~420 cm⁻¹ |
| H₂O | Bend (~1580 cm⁻¹) | 125 cm⁻¹ |

## Layout

geometries/ # input XYZ files
runs/ # MD trajectories (dipole.dat, energy.dat, …)
results/ # figures (.png)
scripts/ # post-processing (compare_bend_hanning.py)
docs/ # report and meeting update (.tex + .pdf)


## Quick start

```bash
cd runs/h2o
python ../../scripts/compare_bend_hanning.py
```

Produces `h2o_bend_hanning.png` and `h2o_full_hanning.png`.
