#!/bin/bash
# Check that the public cboamd (--code cboamd) reproduces runs made earlier with
# polaritonic_deep_md. Each run copies the settings of a 2026-10-02 run exactly:
#   CO2 NEP, omega 2400, lambda 0.1, 10000 steps, chi on and off
#       vs runs/co2/2026-10-02_lambda_sweep/nep_cav_w2400_lam0.1_{chi,nochi}
#   H2O PySCF, omega 1584, lambda 0.1, 800 steps, chi off and on
#       vs runs/h2o/2026-10-02_lambda_sweep/pyscf_cav_w1584_lam0.1_{nochi,chi}
# All four run in parallel; the H2O chi-on run is the slowest (~11 min).
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/code_check.sh
# Then compare:
#   python scripts/compare_runs.py
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-11-10_code_check}
export OMP_NUM_THREADS=1
python scripts/run.py h2o pyscf --cavity --lam 0.1 --steps 800 --campaign $CAMPAIGN --name chi --set polar=true &
python scripts/run.py h2o pyscf --cavity --lam 0.1 --steps 800 --campaign $CAMPAIGN --name nochi &
python scripts/run.py co2 nep --cavity --omega 2400 --lam 0.1 --steps 10000 --campaign $CAMPAIGN --name chi &
python scripts/run.py co2 nep --cavity --omega 2400 --lam 0.1 --steps 10000 --campaign $CAMPAIGN --name nochi --set polar=false &
wait
