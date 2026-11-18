#!/bin/bash
# Two-water separation scan: from which separation R do two waters couple only
# through the cavity? (bare bend unchanged, Rabi splitting = sqrt(2) x one water)
#
#   2 H2O, aligned (C2 axes along x), head-to-tail along x, R = 4, 6, 8, 10 A
#   geometries/h2o_N2_headtail-x_R<R>A.xyz  (scripts/make_cluster.py h2o 2 4 6 8 10)
#   PySCF LDA/cc-pVDZ, chi neglected, 800 steps (83 cm^-1 resolution), photon starts at q = 0
#   cavity: bend, 1584 cm^-1, polarized along x, lambda = 0.2
#           single-water reference: runs/h2o/2026-10-02_lambda_sweep, Rabi 239 cm^-1,
#           so sqrt(2) x = 338 cm^-1
#
# Cost (8-core laptop, OMP_NUM_THREADS=1): ~1.2 s/step per run alone -> ~17 min for
# 800 steps; four in parallel ~17-22 min. Run as two batches of four:
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_N2_separation_scan.sh cavity   # 4 cavity runs, ~20 min
#   bash scripts/campaigns/h2o_N2_separation_scan.sh bare     # 4 no-cavity runs, ~20 min
#   bash scripts/campaigns/h2o_N2_separation_scan.sh all      # both batches, ~40 min
# Then:
#   python scripts/plot_h2o_N2_separation.py runs/h2o_N2/2026-11-18_separation_scan
#   python scripts/molecule_motion.py runs/h2o_N2/2026-11-18_separation_scan
#   python scripts/check_energy.py runs/h2o_N2/2026-11-18_separation_scan
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-11-18_separation_scan}
LAM=${LAM:-0.2}
SEPARATIONS="4.0 6.0 8.0 10.0"
RUN="python scripts/run.py h2o_N2 pyscf --steps 800 --campaign $CAMPAIGN"
export OMP_NUM_THREADS=1

WHAT=${1:-all}
for R in $SEPARATIONS; do
    [ -f "geometries/h2o_N2_headtail-x_R${R}A.xyz" ] || python scripts/make_cluster.py h2o 2 "$R"
done
if [[ $WHAT == cavity || $WHAT == all ]]; then
    for R in $SEPARATIONS; do
        $RUN --cavity --lam "$LAM" --set "xyz_file=geometries/h2o_N2_headtail-x_R${R}A.xyz" --name "R${R}" &
    done
    wait
fi
if [[ $WHAT == bare || $WHAT == all ]]; then
    for R in $SEPARATIONS; do
        $RUN --set "xyz_file=geometries/h2o_N2_headtail-x_R${R}A.xyz" --name "R${R}_bare" &
    done
    wait
fi
