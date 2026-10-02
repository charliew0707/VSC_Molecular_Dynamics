#!/bin/bash
# H2O lambda sweep: one cavity mode resonant with the bend (omega = 1584 cm^-1),
# polarized along x. The geometry (geometries/h2o-single-x.xyz) has the C2 axis along
# x, so the bend's dipole oscillates along x, the component `infrared -m ase` uses.
#
#   PySCF : LDA (lda,vwn) / cc-pVDZ, chi neglected, 800 steps (~83 cm^-1 resolution)
#           ~0.17 s/step -> all runs in parallel take ~5 min
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_lambda_sweep.sh
# Then plot:
#   python scripts/plot_h2o_sweep.py runs/h2o/2026-10-02_lambda_sweep
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-10-02_lambda_sweep}
RUN="python scripts/run.py h2o pyscf"

export OMP_NUM_THREADS=1
$RUN --steps 800 --campaign $CAMPAIGN --name bare &
for lam in 0.05 0.1 0.15 0.2 0.25 0.3; do
    $RUN --cavity --lam $lam --steps 800 --campaign $CAMPAIGN --name nochi &
done
wait
