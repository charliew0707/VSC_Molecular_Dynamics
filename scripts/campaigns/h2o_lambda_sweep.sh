#!/bin/bash
# H2O lambda sweep: one cavity mode resonant with the bend (omega = 1584 cm^-1),
# polarized along x. The geometry (geometries/h2o-single-x.xyz) has the C2 axis along
# x, so the bend's dipole oscillates along x, the component `infrared -m ase` uses.
#
#   PySCF : LDA (lda,vwn) / cc-pVDZ, 800 steps (~83 cm^-1 resolution)
#     nochi : chi neglected, ~0.17 s/step -> all runs in parallel ~7 min
#     chi   : chi included (polarizability every step), ~0.8 s/step -> run in two
#             batches of <= 3 so each stays under ~15 min
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_lambda_sweep.sh [nochi|chi1|chi2|all]
# Then plot:
#   python scripts/plot_h2o_sweep.py runs/h2o/2026-10-02_lambda_sweep
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-10-02_lambda_sweep}
RUN="python scripts/run.py h2o pyscf"

WHAT=${1:-all}
export OMP_NUM_THREADS=1
if [[ $WHAT == nochi || $WHAT == all ]]; then
    $RUN --steps 800 --campaign $CAMPAIGN --name bare &
    for lam in 0.05 0.1 0.15 0.2 0.25 0.3; do
        $RUN --cavity --lam $lam --steps 800 --campaign $CAMPAIGN --name nochi &
    done
    wait
fi
chi_batch() {
    for lam in "$@"; do
        $RUN --cavity --lam $lam --steps 800 --campaign $CAMPAIGN --name chi --set polar=true &
    done
    wait
}
if [[ $WHAT == chi1 || $WHAT == all ]]; then chi_batch 0.1 0.2 0.3; fi
if [[ $WHAT == chi2 || $WHAT == all ]]; then chi_batch 0.15 0.25; fi
