#!/bin/bash
# Relaxed photon start: does the q = 0 start explain the extra two-water splitting?
#
# Every earlier run starts the photon at q = 0, where the cavity field
# Ea = omega q - lambda mu is not zero, so the cavity is "plucked" at t = 0 with
# energy ~ (lambda mu)^2 / 2, four times larger for two waters than for one.
# Here the photon starts relaxed, q0 = lambda mu_x / omega (Ea = 0 at t = 0), with
# mu_x the starting dipole (from the q = 0 runs' dipole.dat):
#   1 H2O: mu_x = -0.774374 a.u. -> q0 = -21.4590
#   2 H2O (R = 8 A): mu_x = -1.553619 a.u. -> q0 = -43.0530
# omega = 1584 cm^-1 = 0.00721724 a.u., lambda = 0.2, chi off, 800 steps.
# Compare with the q = 0 runs:
#   runs/h2o/2026-10-02_lambda_sweep/pyscf_cav_w1584_lam0.2_nochi       (Rabi 239)
#   runs/h2o_N2/2026-11-18_separation_scan/pyscf_cav_w1584_lam0.2_R8.0  (Rabi 375)
# Cost: ~4 min (1 H2O) and ~17 min (2 H2O), run in parallel.
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_relaxed_photon.sh
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-11-18_relaxed_photon}
export OMP_NUM_THREADS=1
python scripts/run.py h2o pyscf --cavity --lam 0.2 --steps 800 --campaign $CAMPAIGN \
    --set qa0=-21.4590 --name relaxed &
python scripts/run.py h2o_N2 pyscf --cavity --lam 0.2 --steps 800 --campaign $CAMPAIGN \
    --set xyz_file=geometries/h2o_N2_headtail-x_R8.0A.xyz --set qa0=-43.0530 --name R8.0_relaxed &
wait
