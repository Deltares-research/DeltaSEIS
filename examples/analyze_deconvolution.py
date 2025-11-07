"""
Analyze the deconvolution filter response to understand the notch issue
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Simulate what happens in the deconvolution

# Use the extracted wavelet properties from previous analysis
# Trace 6415 (0-indexed): Peak at 5068 Hz
# Let's create a synthetic wavelet spectrum similar to what we observed

fs = 50_000  # Hz
n_samples = 3256
dt = 1/fs
freqs = np.fft.rfftfreq(n_samples, dt)

# Create a simplified model of the signature spectrum
# Peaked around 5000 Hz with some bandwidth
center_freq = 5068
bandwidth = 1000
sig_amplitude = np.exp(-((freqs - center_freq)**2) / (2 * (bandwidth/2.355)**2))
sig_amplitude = sig_amplitude / np.max(sig_amplitude)  # normalize

# Add some realistic noise floor
sig_amplitude = sig_amplitude * 0.95 + 0.05

# Convert to complex FFT (assuming minimum phase for simplicity)
sig_fft = sig_amplitude.astype(complex)
sig_power = np.abs(sig_fft)**2

# Test different Wiener filter formulations
epsilon_values = [0.001, 0.01, 0.05, 0.1]

fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle('Wiener Filter Analysis - Different Stabilization Methods', fontsize=14, fontweight='bold')

# Method 1: Current implementation (proportional stabilization)
ax = axes[0, 0]
for eps in epsilon_values:
    stabilization = eps * sig_power
    wiener_filter = np.conj(sig_fft) / (sig_power + stabilization)
    gain = np.abs(wiener_filter) * np.abs(sig_fft)  # Effective gain
    ax.plot(freqs/1000, gain, label=f'ε={eps}', linewidth=1.5)
ax.set_title('Method 1: ε × sig_power (Current)')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Effective Gain')
ax.set_xlim([0, 15])
ax.legend()
ax.grid(True, alpha=0.3)
ax.axvline(center_freq/1000, color='r', linestyle='--', alpha=0.3, label='Peak freq')

# Method 2: Constant stabilization
ax = axes[0, 1]
for eps in epsilon_values:
    stabilization = eps * np.mean(sig_power)  # Constant
    wiener_filter = np.conj(sig_fft) / (sig_power + stabilization)
    gain = np.abs(wiener_filter) * np.abs(sig_fft)
    ax.plot(freqs/1000, gain, label=f'ε={eps}', linewidth=1.5)
ax.set_title('Method 2: ε × mean(sig_power)')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Effective Gain')
ax.set_xlim([0, 15])
ax.legend()
ax.grid(True, alpha=0.3)
ax.axvline(center_freq/1000, color='r', linestyle='--', alpha=0.3)

# Method 3: Power-dependent but with exponent
ax = axes[1, 0]
for eps in epsilon_values:
    stabilization = eps * np.sqrt(sig_power) * np.mean(np.sqrt(sig_power))
    wiener_filter = np.conj(sig_fft) / (sig_power + stabilization)
    gain = np.abs(wiener_filter) * np.abs(sig_fft)
    ax.plot(freqs/1000, gain, label=f'ε={eps}', linewidth=1.5)
ax.set_title('Method 3: ε × √sig_power × mean(√sig_power)')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Effective Gain')
ax.set_xlim([0, 15])
ax.legend()
ax.grid(True, alpha=0.3)
ax.axvline(center_freq/1000, color='r', linestyle='--', alpha=0.3)

# Method 4: Adaptive (stronger at high power)
ax = axes[1, 1]
for eps in epsilon_values:
    # Use geometric mean for stability
    stabilization = eps * np.exp(np.mean(np.log(sig_power + 1e-10)))
    wiener_filter = np.conj(sig_fft) / (sig_power + stabilization)
    gain = np.abs(wiener_filter) * np.abs(sig_fft)
    ax.plot(freqs/1000, gain, label=f'ε={eps}', linewidth=1.5)
ax.set_title('Method 4: ε × exp(mean(log(sig_power)))')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Effective Gain')
ax.set_xlim([0, 15])
ax.legend()
ax.grid(True, alpha=0.3)
ax.axvline(center_freq/1000, color='r', linestyle='--', alpha=0.3)

# Show the signature power spectrum
ax = axes[2, 0]
ax.semilogy(freqs/1000, sig_power, 'b-', linewidth=2, label='Signature power')
ax.axvline(center_freq/1000, color='r', linestyle='--', label=f'Peak: {center_freq} Hz')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Power (log scale)')
ax.set_title('Input Signature Power Spectrum')
ax.set_xlim([0, 15])
ax.legend()
ax.grid(True, alpha=0.3)

# Show ratio sig_power / (sig_power + stabilization) for Method 1
ax = axes[2, 1]
for eps in epsilon_values:
    stabilization = eps * sig_power
    ratio = sig_power / (sig_power + stabilization)
    ax.plot(freqs/1000, ratio, label=f'ε={eps}', linewidth=1.5)
ax.set_title('Method 1: Attenuation Factor\n(lower = more suppression)')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('sig_power / (sig_power + ε×sig_power)')
ax.set_xlim([0, 15])
ax.set_ylim([0, 1.05])
ax.legend()
ax.grid(True, alpha=0.3)
ax.axvline(center_freq/1000, color='r', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('deconvolution_filter_analysis.png', dpi=150, bbox_inches='tight')
print("Saved: deconvolution_filter_analysis.png")
plt.show()

print("\n=== KEY INSIGHTS ===")
print(f"Signature peak frequency: {center_freq} Hz")
print(f"\nFor Method 1 (current: ε × sig_power):")
print("  - Attenuation factor = sig_power / (sig_power + ε×sig_power) = 1/(1+ε)")
print("  - This is CONSTANT across all frequencies!")
print("  - At peak: suppressed by factor 1/(1+ε)")
print("  - At weak frequencies: also suppressed by factor 1/(1+ε)")
print("  - Problem: NO frequency-dependent stabilization - not what we want!")
print(f"\nWith ε=0.01: gain is {1/(1+0.01):.4f} everywhere (99% of signal)")
print(f"With ε=0.1: gain is {1/(1+0.1):.4f} everywhere (91% of signal)")

print("\n=== RECOMMENDATION ===")
print("Method 1 is theoretically correct for Wiener filtering!")
print("The 'notch' you're seeing might be due to:")
print("  1. Phase issues (not addressed by amplitude-only Wiener)")
print("  2. The signature itself having strong power at 5kHz")
print("  3. Need for different epsilon value")
print("  4. Need for additional bandpass filtering post-decon")
