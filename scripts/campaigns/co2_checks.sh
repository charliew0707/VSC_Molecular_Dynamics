#!/bin/bash
# CO2 single-molecule checks before multi-molecule runs (NEP, chi off, 10000 steps
# = 6.7 cm^-1 resolution, lambda = 0.1). CO2 has no permanent dipole, so unlike
# H2O the tilted-polarization runs should not torque the molecule.
#
#   orientation : polarization tilted from x (the CO2 axis) toward z by theta.
#                 Expect Omega_R ~ cos(theta).
#   detuning    : polarization along x, cavity frequency swept across the
#                 asymmetric stretch (bare NEP stretch: 2438 cm^-1).
#   amplitude   : initial C displacement along x (asymmetric stretch) 0.005-0.04 A;
#                 geometries/co2-disp.xyz is the 0.01 A case.
#
# Usage (from the repo root, in the vsc env):  bash scripts/campaigns/co2_checks.sh
# Then: python scripts/plot_co2_checks.py
set -e
cd "$(dirname "$0")/../.."
RUN="python scripts/run.py co2 nep --cavity --lam 0.1 --steps 10000 --set polar=false"
export OMP_NUM_THREADS=1

for th in 0 15 30 45 60 75 90; do
    vec=$(python -c "import math; t=math.radians($th); print(f'[[{math.cos(t):.6f}, 0, {math.sin(t):.6f}]]')")
    $RUN --omega 2400 --campaign 2026-10-02_orientation --name th$th --set "lambda_vector=$vec" &
done
for w in 2100 2200 2300 2400 2438 2500 2600 2700 2800; do
    $RUN --omega $w --campaign 2026-10-02_detuning --name detune &
done
wait
for a in 0.005 0.01 0.02 0.04; do
    geom=geometries/co2-disp-$a.xyz; [[ $a == 0.01 ]] && geom=geometries/co2-disp.xyz
    $RUN --omega 2400 --campaign 2026-10-02_amplitude --name disp$a --set xyz_file=$geom &
done
wait
