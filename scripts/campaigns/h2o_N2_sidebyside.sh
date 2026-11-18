#!/bin/bash
# Two waters side by side: the head-to-tail pair (R = 8 A) attracts and collides at
# ~1.35 ps (scripts/campaigns/h2o_relaxed_photon_long.sh). Side by side, the two
# parallel dipoles repel, so the pair should stay apart for the whole run.
#
#   2 H2O, aligned (C2 axes along x), shifted along y by R = 8 A
#   geometries/h2o_N2_sidebyside-y_R8.0A.xyz (scripts/make_cluster.py h2o 2 8.0 --arrangement sidebyside)
#   PySCF LDA/cc-pVDZ, chi off, cavity at the bend (1584 cm^-1) polarized along x, lambda 0.2
#   relaxed photon start q0 = lambda mu_x / omega, starting dipole from a PySCF single point:
#     mu_x = -1.545982 a.u. (head-to-tail: -1.553629; one water: -0.774379)
#     omega = 0.00721724 a.u. -> q0 = -42.8414
#   Compare with 1 H2O relaxed, 4000 steps: Rabi 250 cm^-1, so sqrt(2) x = 354.
#
# Results, 800 steps (2026-11-18, cboamd 89bc6dd):
#   bare:   bend 1582 cm^-1 (same as head-to-tail R = 8 A and one water)
#   cavity: LP 1426 / UP 1770, Rabi 344 (head-to-tail relaxed, 800 steps: 354)
#   motion: the pair pushes apart (separation 8.00 -> 8.21 A) and each water turns
#           by 15 deg, with or without the cavity (dipole-dipole torque), so
#           check the rotation in the long run.
#   energy ok in both (fluct 1.9%, drift < 0.3%).
# Results, 4000 steps (2 ps, 8 cm^-1 resolution):
#   cavity: LP 1424 / UP 1764, Rabi 340 (1st ps 334, 2nd ps 350); ratio to one
#           water (250) = 1.36, vs sqrt(2) = 1.41. The waters librate (each swings
#           up to 22 deg from x and back, mirror images of each other); average
#           cos ~ 0.97 accounts for most of the 4% shortfall.
#   motion: the pair keeps pushing apart, 8.0 -> 11.6 A; no contact.
#   energy: check_energy says CHECK (scale only 16.5 meV), but the absolute H
#           range is 1.0 meV and the drift -0.5 meV over 2 ps: fine.
#
# Usage (from the repo root, in the vsc env):
#   bash scripts/campaigns/h2o_N2_sidebyside.sh short   # cavity + bare, 800 steps, ~18 min
#   bash scripts/campaigns/h2o_N2_sidebyside.sh long    # cavity, 4000 steps, ~1.5 h (own terminal)
# Then:
#   python scripts/molecule_motion.py runs/h2o_N2/2026-11-18_sidebyside*
#   python scripts/check_energy.py runs/h2o_N2/2026-11-18_sidebyside*
set -e
cd "$(dirname "$0")/../.."
XYZ=geometries/h2o_N2_sidebyside-y_R8.0A.xyz
Q0=-42.8414
export OMP_NUM_THREADS=1
[ -f $XYZ ] || python scripts/make_cluster.py h2o 2 8.0 --arrangement sidebyside
case ${1:-short} in
short)
    python scripts/run.py h2o_N2 pyscf --cavity --lam 0.2 --steps 800 --campaign 2026-11-18_sidebyside \
        --set xyz_file=$XYZ --set qa0=$Q0 --name sidebyside_R8.0_relaxed &
    python scripts/run.py h2o_N2 pyscf --steps 800 --campaign 2026-11-18_sidebyside \
        --set xyz_file=$XYZ --name sidebyside_R8.0_bare &
    wait ;;
long)
    python scripts/run.py h2o_N2 pyscf --cavity --lam 0.2 --steps 4000 --campaign 2026-11-18_sidebyside_long \
        --set xyz_file=$XYZ --set qa0=$Q0 --name sidebyside_R8.0_relaxed ;;
esac
