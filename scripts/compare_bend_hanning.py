"""
compare_bend_hanning.py
Run from: runs/h2o/  (reads pyscf/dipole.dat and pyscf_cavity/dipole.dat)

Method:
- Reads raw dipole.dat instead of spectrum_ase.dat
- Applies a Hanning window before FFT to kill spectral leakage
- Zero-pads 4x for a smoother frequency axis
- Shows stacked no-cavity / with-cavity plots with peak annotations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

TIMESTEP_FS = 0.5          # must match in.json timestep
C_CM_S      = 2.99792458e10  # speed of light in cm/s

def dipole_to_spectrum(dipole_dat_path, col=7, fmin=800, fmax=4500):
    """
    col 5/6/7 = total dipole x/y/z  (what you want for IR)
    col 2/3/4 = dressed dipole x/y/z (includes cavity field; only in cavity runs)
    """
    data = np.loadtxt(dipole_dat_path)
    mu   = data[:, col]
    N    = len(mu)

    # Hanning window — tapers the signal to zero at both ends,
    # which removes the 'ringing' artefacts that appear when you
    # just cut off a trajectory abruptly.
    window = np.hanning(N)
    mu_w   = mu * window

    # Zero-pad 4x → finer frequency grid (does NOT add resolution,
    # but makes the curve look smoother — think of it as interpolation)
    N_pad  = 4 * N
    fft    = np.fft.rfft(mu_w, n=N_pad)
    power  = np.abs(fft) ** 2

    dt_s   = TIMESTEP_FS * 1e-15          # timestep in seconds
    freq_hz = np.fft.rfftfreq(N_pad, d=dt_s)
    freq   = freq_hz / C_CM_S             # convert to cm⁻¹

    mask   = (freq > fmin) & (freq < fmax)
    return freq[mask], power[mask]


def normalise(arr, fmin, fmax, freq):
    """Normalise to the max within a sub-window [fmin, fmax]."""
    sub = arr[(freq > fmin) & (freq < fmax)]
    mx  = sub.max() if sub.max() > 0 else 1.0
    return arr / mx


# ── Load ──────────────────────────────────────────────────────────────────────
freq_bare, pow_bare = dipole_to_spectrum("pyscf/dipole.dat",         col=7)
freq_cav,  pow_cav  = dipole_to_spectrum("pyscf_cavity/dipole.dat",  col=7)

# Normalise each spectrum to its max in the bend region (1200–2000 cm⁻¹)
pow_bare_n = normalise(pow_bare, 1200, 2000, freq_bare)
pow_cav_n  = normalise(pow_cav,  1200, 2000, freq_cav)

# ── Plot 1: bend region stacked ───────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
fig.suptitle("H₂O bend mode: Rabi splitting with Hanning window", fontsize=13)

# Mask to bend region for display
b_mask = (freq_bare > 1200) & (freq_bare < 2000)
c_mask = (freq_cav  > 1200) & (freq_cav  < 2000)

axes[0].plot(freq_bare[b_mask], pow_bare_n[b_mask], color='dimgray', lw=1.5)
axes[0].axvline(1595, color='green', lw=1, ls=':', alpha=0.8, label='Expt. bend 1595 cm⁻¹')
axes[0].set_ylabel("IR intensity (norm.)")
axes[0].set_title("No cavity")
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.25)

axes[1].plot(freq_cav[c_mask], pow_cav_n[c_mask], color='steelblue', lw=1.5)
axes[1].axvline(1584, color='orange', lw=1, ls=':', alpha=0.8, label='Cavity ω = 1584 cm⁻¹')

peaks, _ = find_peaks(pow_cav_n[c_mask], height=0.1, distance=8, prominence=0.06)
fc = freq_cav[c_mask]
pc = pow_cav_n[c_mask]
for p in peaks:
    axes[1].axvline(fc[p], color='red', lw=0.9, ls=':', alpha=0.85)
    axes[1].annotate(f"{fc[p]:.0f}", xy=(fc[p], pc[p]),
                     xytext=(fc[p]+18, pc[p]-0.12), fontsize=9, color='red',
                     arrowprops=dict(arrowstyle='->', color='red', lw=0.7))

axes[1].set_xlabel("Frequency (cm⁻¹)")
axes[1].set_ylabel("IR intensity (norm.)")
axes[1].set_title("With cavity  (λ = 0.1, z-polarisation)")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig("h2o_bend_hanning.png", dpi=150)
plt.close()
print("✓  Saved  h2o_bend_hanning.png")

# ── Plot 2: full spectrum (all modes) ─────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.plot(freq_bare, pow_bare / pow_bare.max(), color='dimgray',
         lw=1.2, label='No cavity', alpha=0.85)
ax2.plot(freq_cav,  pow_cav  / pow_cav.max(),  color='steelblue',
         lw=1.2, label='Cavity (λ=0.1)', alpha=0.85)
for x, lbl in [(1595, 'bend'), (3657, 'sym-str'), (3756, 'asym-str')]:
    ax2.axvline(x, color='green', lw=0.8, ls='--', alpha=0.5)
    ax2.text(x+25, 0.92, lbl, color='green', fontsize=8)
ax2.set_xlabel("Frequency (cm⁻¹)")
ax2.set_ylabel("IR intensity (norm. globally)")
ax2.set_title("H₂O full spectrum — no cavity vs cavity")
ax2.legend()
ax2.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig("h2o_full_hanning.png", dpi=150)
plt.close()
print("✓  Saved  h2o_full_hanning.png")

# ── Console summary ───────────────────────────────────────────────────────────
print(f"\nTrajectory length: {len(np.loadtxt('pyscf/dipole.dat'))} steps  →  "
      f"Δν ≈ {66667/len(np.loadtxt('pyscf/dipole.dat')):.0f} cm⁻¹ raw resolution")

if len(peaks) >= 2:
    lp = fc[peaks[0]]
    up = fc[peaks[-1]]
    mid = 0.5 * (lp + up)
    split = up - lp
    print(f"\nCavity peaks found:")
    print(f"  LP  =  {lp:.0f} cm⁻¹")
    print(f"  UP  =  {up:.0f} cm⁻¹")
    print(f"  mid =  {mid:.0f} cm⁻¹  (cavity at 1584)")
    print(f"  Rabi splitting ≈ {split:.0f} cm⁻¹")
elif len(peaks) == 1:
    print(f"\nOnly one cavity peak found at {fc[peaks[0]]:.0f} cm⁻¹ — "
          "try running longer (e.g. steps=5000) for cleaner splitting.")
else:
    print("\nNo peaks detected above threshold — spectrum may need more steps.")
