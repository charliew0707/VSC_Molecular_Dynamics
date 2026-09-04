import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def load_spec(path, fmin=1200, fmax=2000):
    d = np.loadtxt(path)
    freq, intensity = d[:,0], d[:,1]
    mask = (freq > fmin) & (freq < fmax)
    freq_r = freq[mask]
    intensity_r = intensity[mask]
    return freq_r, intensity_r / intensity_r.max()

freq_bare, spec_bare = load_spec("pyscf/spectrum_ase.dat")
freq_cav,  spec_cav  = load_spec("pyscf_cavity/spectrum_ase.dat")

fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(freq_bare, spec_bare, color='gray',      lw=1.4, linestyle='--', label='No cavity')
ax.plot(freq_cav,  spec_cav,  color='steelblue', lw=1.4, label='Cavity (λ=0.1, ω=1584 cm⁻¹)')
ax.axvline(1584, color='orange', lw=1, linestyle=':', alpha=0.7, label='Cavity freq (1584 cm⁻¹)')

peaks, _ = find_peaks(spec_cav, height=0.15, distance=10, prominence=0.1)
for p in peaks:
    ax.axvline(freq_cav[p], color='red', lw=0.8, linestyle=':', alpha=0.8)
    ax.annotate(f"{freq_cav[p]:.0f}", xy=(freq_cav[p], spec_cav[p]),
                xytext=(freq_cav[p]+15, spec_cav[p]-0.1), fontsize=9, color='red')

ax.set_xlabel("Frequency (cm⁻¹)")
ax.set_ylabel("IR Intensity (norm. to local max)")
ax.set_title("H₂O bend mode — cavity vs no cavity (PySCF)")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("h2o_bend_comparison.png", dpi=150)
print("Saved h2o_bend_comparison.png")
if len(peaks) >= 2:
    print(f"LP: {freq_cav[peaks[0]]:.0f} cm⁻¹, UP: {freq_cav[peaks[-1]]:.0f} cm⁻¹")
    print(f"Rabi splitting: {freq_cav[peaks[-1]] - freq_cav[peaks[0]]:.0f} cm⁻¹")
else:
    print(f"Peaks found: {[f'{freq_cav[p]:.0f}' for p in peaks]} cm⁻¹")
