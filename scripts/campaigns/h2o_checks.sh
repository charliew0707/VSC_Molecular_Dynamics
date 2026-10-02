#!/bin/bash
# H2O single-molecule checks before multi-molecule runs (PySCF LDA/cc-pVDZ, chi off,
# 800 steps, lambda = 0.3 so splittings stay resolvable at 83 cm^-1).
#
#   orientation : cavity polarization tilted from x (the C2 axis / bend dipole)
#                 toward z (out of the molecular plane) by theta; the molecule is
#                 not rotated, so the bend dipole stays along x for `infrared`.
#                 Expect Omega_R ~ cos(theta).
#   detuning    : polarization along x, cavity frequency swept across the bend.
#                 Expect an anti-crossing.
#
# Usage (from the repo root, in the vsc env):  bash scripts/campaigns/h2o_checks.sh
# Then: python scripts/plot_h2o_checks.py
set -e
cd "$(dirname "$0")/../.."
LAM=0.3
RUN="python scripts/run.py h2o pyscf --cavity --lam $LAM --steps 800"
export OMP_NUM_THREADS=1

for th in 0 30 45 60 75 90; do
    vec=$(python -c "import math; t=math.radians($th); print(f'[[{math.cos(t):.6f}, 0, {math.sin(t):.6f}]]')")
    $RUN --campaign 2026-10-02_orientation --name th$th --set "lambda_vector=$vec" &
done
for w in 1300 1400 1500 1584 1650 1750 1850; do
    $RUN --omega $w --campaign 2026-10-02_detuning --name detune &
done
wait
