#!/bin/bash
# CO2 lambda sweep, benchmark-style (cf. polaritonic_deep_md/tests/e2e/compare_with_bonini2024):
# one cavity mode at omega = 2400 cm^-1 polarized along x (the CO2 axis), lambda swept.
#
#   NEP   : 10000 steps (6.7 cm^-1 resolution), polarizability (chi) included AND neglected
#   PySCF : 800 steps (84 cm^-1 resolution; keeps the parallel batch under ~25 min), chi neglected (chi on costs ~2x per step)
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/co2_lambda_sweep.sh [nep|pyscf|all]
# Then plot:
#   python scripts/plot_co2_sweep.py runs/co2/2026-10-02_lambda_sweep
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-10-02_lambda_sweep}
WHAT=${1:-all}
RUN="python scripts/run.py co2"

if [[ $WHAT == nep || $WHAT == all ]]; then
    LAMS="0.02 0.05 0.1 0.13 0.17 0.2 0.25 0.3"   # subset of the benchmark's lambda grid
    $RUN nep --steps 10000 --campaign $CAMPAIGN --name bare
    $RUN nep --cavity --omega 2400 --lam $LAMS --steps 10000 --campaign $CAMPAIGN --name chi
    $RUN nep --cavity --omega 2400 --lam $LAMS --steps 10000 --campaign $CAMPAIGN --name nochi --set polar=false
fi

if [[ $WHAT == pyscf || $WHAT == all ]]; then
    # one single-threaded process per lambda, all in parallel (~23 min on 8 cores)
    export OMP_NUM_THREADS=1
    $RUN pyscf --steps 800 --campaign $CAMPAIGN --name bare &
    for lam in 0.05 0.1 0.2 0.3; do
        $RUN pyscf --cavity --omega 2400 --lam $lam --steps 800 --campaign $CAMPAIGN --name nochi &
    done
    wait
fi
