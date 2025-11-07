"""
Test deconvolution with different parameters on actual data
"""
import numpy as np
import matplotlib.pyplot as plt
import deltaseis
from pathlib import Path

# Load the data
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]
segy_path = segy_paths[0]

print(f"Loading: {segy_path.name}")
segy = deltaseis.Segy_edit(segy_path)
segy.xstar_split('full')

trace_array = np.array(segy.trace_data)
if trace_array.shape[0] == 12404:
    trace_array = trace_array.T

# Test on a smaller subset for speed
trace_subset = trace_array[:, 6000:6800]  # 800 traces around the signature trace

fs = 50_000
dt = 1/fs

# Create comparison plot
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Deconvolution Parameter Comparison', fontsize=14, fontweight='bold')

# Original data
seis_orig = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_orig.time_power_gain(3)
seis_orig.trace_averaging(1)

# Plot original
ax = axes[0, 0]
time_ms = np.arange(trace_subset.shape[0]) * dt * 1000
extent = [0, trace_subset.shape[1], time_ms[-1], time_ms[0]]
vmax = np.percentile(np.abs(seis_orig.data), 98)
im = ax.imshow(seis_orig.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('Original (no decon)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([20, 45])  # Focus on area around 30-32 ms

# Compute and plot frequency spectrum of original
fft_orig = np.fft.rfft(seis_orig.data, axis=0)
freqs = np.fft.rfftfreq(seis_orig.data.shape[0], dt)
power_orig = np.mean(np.abs(fft_orig)**2, axis=1)
ax_spec = axes[1, 0]
ax_spec.semilogy(freqs/1000, power_orig, 'b-', linewidth=1.5, label='Original')
ax_spec.set_xlabel('Frequency (kHz)')
ax_spec.set_ylabel('Average Power')
ax_spec.set_xlim([0, 20])
ax_spec.grid(True, alpha=0.3)
ax_spec.legend()
ax_spec.set_title('Frequency Content')

# Test 1: Current method with epsilon=0.01
seis1 = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis1.time_power_gain(3)
seis1.trace_averaging(1)
# Note: trace 6415 in full data = trace 415 in subset (6415 - 6000)
seis1.signature_deconvolution(415, 30.5, 31.7, method='wiener', epsilon=0.01, prewhiten=True, prewhiten_percent=1.0)

ax = axes[0, 1]
im = ax.imshow(seis1.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('Decon: ε=0.01, prewhiten=1%')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([20, 45])

fft1 = np.fft.rfft(seis1.data, axis=0)
power1 = np.mean(np.abs(fft1)**2, axis=1)
ax_spec = axes[1, 1]
ax_spec.semilogy(freqs/1000, power_orig, 'b-', linewidth=1, alpha=0.5, label='Original')
ax_spec.semilogy(freqs/1000, power1, 'r-', linewidth=1.5, label='Decon')
ax_spec.axvline(5.068, color='orange', linestyle='--', alpha=0.7, label='Wavelet peak')
ax_spec.set_xlabel('Frequency (kHz)')
ax_spec.set_ylabel('Average Power')
ax_spec.set_xlim([0, 20])
ax_spec.grid(True, alpha=0.3)
ax_spec.legend()
ax_spec.set_title('Frequency - check for notch at 5kHz')

# Test 2: Larger epsilon
seis2 = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis2.time_power_gain(3)
seis2.trace_averaging(1)
seis2.signature_deconvolution(415, 30.5, 31.7, method='wiener', epsilon=0.1, prewhiten=True, prewhiten_percent=1.0)

ax = axes[0, 2]
im = ax.imshow(seis2.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('Decon: ε=0.1 (more conservative)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([20, 45])

fft2 = np.fft.rfft(seis2.data, axis=0)
power2 = np.mean(np.abs(fft2)**2, axis=1)
ax_spec = axes[1, 2]
ax_spec.semilogy(freqs/1000, power_orig, 'b-', linewidth=1, alpha=0.5, label='Original')
ax_spec.semilogy(freqs/1000, power2, 'g-', linewidth=1.5, label='Decon')
ax_spec.axvline(5.068, color='orange', linestyle='--', alpha=0.7, label='Wavelet peak')
ax_spec.set_xlabel('Frequency (kHz)')
ax_spec.set_ylabel('Average Power')
ax_spec.set_xlim([0, 20])
ax_spec.grid(True, alpha=0.3)
ax_spec.legend()
ax_spec.set_title('Frequency - check for notch at 5kHz')

plt.tight_layout()
plt.savefig('deconvolution_parameter_test.png', dpi=150, bbox_inches='tight')
print("\nSaved: deconvolution_parameter_test.png")
plt.show()

print("\n=== NOTCH DETECTION ===")
# Check if there's a notch at 5kHz
idx_5khz = np.argmin(np.abs(freqs - 5068))
freq_window = slice(idx_5khz-10, idx_5khz+10)

power_orig_5k = np.mean(power_orig[freq_window])
power1_5k = np.mean(power1[freq_window])
power2_5k = np.mean(power2[freq_window])

print(f"Average power around 5kHz:")
print(f"  Original: {power_orig_5k:.2e}")
print(f"  Decon (ε=0.01): {power1_5k:.2e} (ratio: {power1_5k/power_orig_5k:.2f}x)")
print(f"  Decon (ε=0.1): {power2_5k:.2e} (ratio: {power2_5k/power_orig_5k:.2f}x)")

if power1_5k < power_orig_5k * 0.5:
    print("\n⚠️ NOTCH DETECTED with ε=0.01: Power suppressed by >50% at 5kHz")
else:
    print("\n✓ No strong notch with ε=0.01")

if power2_5k < power_orig_5k * 0.5:
    print("⚠️ NOTCH DETECTED with ε=0.1: Power suppressed by >50% at 5kHz")
else:
    print("✓ No strong notch with ε=0.1")
