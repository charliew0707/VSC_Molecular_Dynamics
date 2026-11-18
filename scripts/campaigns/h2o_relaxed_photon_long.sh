#!/bin/bash
# Long (4000 steps = 2 ps, 8 cm^-1 resolution) version of h2o_relaxed_photon.sh:
# 1 H2O and 2 H2O (head-to-tail, R = 8 A), lambda = 0.2, chi off, relaxed photon start.
#
# Results (2026-11-18, cboamd 89bc6dd):
#   1 H2O: LP 1464 / UP 1714, Rabi 250 (same in both 1 ps halves; full run 1465/1713)
#   2 H2O: the two waters attract (aligned dipoles, head-to-tail), drift together,
#          touch at ~1.35 ps (closest atoms 2.3 A) and bounce apart. Use only the
#          first 1 ps (steps 0-2000, separation still > 5.8 A): LP 1418 / UP 1772,
#          Rabi 354 -> ratio 354/250 = 1.42 (sqrt 2 = 1.414). The full-run spectrum
#          has doubled peaks (1411/1432, 1759/1782) from the collision.
#   Energy: 2 H2O fluct 0.15%, drift 0.03% (ok). 1 H2O is flagged CHECK only because
#          check_energy normalizes by the tiny moving energy (6.5 meV) of a relaxed
#          start; absolute H range 0.35 meV, drift 0.08 meV over 2 ps.
# Cost: 1 H2O ~21 min; 2 H2O ~93 min (Charlie runs this one in his own terminal).
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_relaxed_photon_long.sh 1   # single water
#   bash scripts/campaigns/h2o_relaxed_photon_long.sh 2   # two waters (~1.5 h)
set -e
cd "$(dirname "$0")/../.."
CAMPAIGN=${CAMPAIGN:-2026-11-18_relaxed_photon_long}
export OMP_NUM_THREADS=1
case ${1:-1} in
1) python scripts/run.py h2o pyscf --cavity --lam 0.2 --steps 4000 --campaign $CAMPAIGN \
       --set qa0=-21.4590 --name relaxed ;;
2) python scripts/run.py h2o_N2 pyscf --cavity --lam 0.2 --steps 4000 --campaign $CAMPAIGN \
       --set xyz_file=geometries/h2o_N2_headtail-x_R8.0A.xyz --set qa0=-43.0530 --name R8.0_relaxed ;;
esac
