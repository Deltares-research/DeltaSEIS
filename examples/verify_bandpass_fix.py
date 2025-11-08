"""
Verify bandpass filter is working correctly after fix
"""
import numpy as np
import matplotlib.pyplot as plt
import deltaseis
from pathlib import Path

# Load data
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]
segy_path = segy_paths[0]

print(f"Loading: {segy_path.name}\n")
s = deltaseis.Segy_edit(segy_path)
s.xstar_split('full')

# Work with subset for speed
trace_subset = np.array(s.trace_data).T[:, 5000:5200]
fs = 50_000
dt = 1/fs

print("="*60)
print("BANDPASS FILTER VERIFICATION TEST")
print("="*60)

# Test 1: Original data
seis_orig = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_orig.time_power_gain(3)

# Test 2: With bandpass filter
seis_filt = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_filt.time_power_gain(3)
seis_filt.bandpass_filter(lowcut=2700, highcut=8000)

# Compute spectra
fft_orig = np.fft.rfft(seis_orig.data, axis=0)
fft_filt = np.fft.rfft(seis_filt.data, axis=0)
freqs = np.fft.rfftfreq(seis_orig.data.shape[0], dt)

power_orig = np.mean(np.abs(fft_orig)**2, axis=1)
power_filt = np.mean(np.abs(fft_filt)**2, axis=1)

# Create comprehensive comparison
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Bandpass Filter Verification (2.7-8.0 kHz)', fontsize=14, fontweight='bold')

# 1. Frequency spectrum - linear
ax = axes[0, 0]
ax.plot(freqs/1000, power_orig, 'b-', linewidth=1.5, alpha=0.7, label='Original')
ax.plot(freqs/1000, power_filt, 'r-', linewidth=1.5, label='Filtered')
ax.axvline(2.7, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Passband')
ax.axvline(8.0, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax.axvspan(2.7, 8.0, alpha=0.1, color='green')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Power')
ax.set_xlim([0, 15])
ax.set_title('Power Spectrum (Linear Scale)')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Frequency spectrum - log
ax = axes[0, 1]
ax.semilogy(freqs/1000, power_orig, 'b-', linewidth=1.5, alpha=0.7, label='Original')
ax.semilogy(freqs/1000, power_filt, 'r-', linewidth=1.5, label='Filtered')
ax.axvline(2.7, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax.axvline(8.0, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax.axvspan(2.7, 8.0, alpha=0.1, color='green', label='Passband')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Power (log scale)')
ax.set_xlim([0, 15])
ax.set_title('Power Spectrum (Log Scale)')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Attenuation
ax = axes[0, 2]
attenuation_db = 10 * np.log10(power_filt / (power_orig + 1e-10))
ax.plot(freqs/1000, attenuation_db, 'purple', linewidth=2)
ax.axhline(0, color='k', linestyle='-', linewidth=1, alpha=0.5)
ax.axhline(-3, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='-3 dB')
ax.axvline(2.7, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Cutoff freq')
ax.axvline(8.0, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax.axvspan(2.7, 8.0, alpha=0.1, color='green')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Attenuation (dB)')
ax.set_xlim([0, 15])
ax.set_ylim([-60, 5])
ax.set_title('Filter Response (dB)')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Time domain - original
ax = axes[1, 0]
time_ms = np.arange(seis_orig.data.shape[0]) * dt * 1000
extent = [0, seis_orig.data.shape[1], time_ms[-1], time_ms[0]]
vmax = np.percentile(np.abs(seis_orig.data), 98)
im = ax.imshow(seis_orig.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('Original Data')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

# 5. Time domain - filtered
ax = axes[1, 1]
im = ax.imshow(seis_filt.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('After Bandpass (2.7-8.0 kHz)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

# 6. Single trace comparison
ax = axes[1, 2]
trace_idx = 100
ax.plot(seis_orig.data[:, trace_idx], time_ms, 'b-', linewidth=1, alpha=0.6, label='Original')
ax.plot(seis_filt.data[:, trace_idx], time_ms, 'r-', linewidth=1.5, label='Filtered')
ax.set_xlabel('Amplitude')
ax.set_ylabel('Time (ms)')
ax.set_title(f'Single Trace Comparison (trace {trace_idx})')
ax.set_ylim([50, 20])
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('bandpass_verification.png', dpi=150, bbox_inches='tight')
print("\nSaved: bandpass_verification.png")
plt.show()

# Quantitative verification
print("\n" + "="*60)
print("QUANTITATIVE VERIFICATION")
print("="*60)

# Check attenuation outside passband
lowfreq_mask = freqs < 2000  # Well below lowcut
highfreq_mask = freqs > 10000  # Well above highcut
passband_mask = (freqs >= 2700) & (freqs <= 8000)

attenuation_low = np.mean(power_filt[lowfreq_mask]) / np.mean(power_orig[lowfreq_mask])
attenuation_high = np.mean(power_filt[highfreq_mask]) / np.mean(power_orig[highfreq_mask])
preservation_pass = np.mean(power_filt[passband_mask]) / np.mean(power_orig[passband_mask])

print(f"\nFrequency < 2 kHz (should be attenuated):")
print(f"  Power ratio: {attenuation_low:.4f} ({10*np.log10(attenuation_low):.1f} dB)")

print(f"\nFrequency > 10 kHz (should be attenuated):")
print(f"  Power ratio: {attenuation_high:.4f} ({10*np.log10(attenuation_high):.1f} dB)")

print(f"\nPassband 2.7-8.0 kHz (should be preserved):")
print(f"  Power ratio: {preservation_pass:.4f} ({10*np.log10(preservation_pass):.1f} dB)")

# Check if data actually changed
data_changed = not np.allclose(seis_orig.data, seis_filt.data)
print(f"\nData was modified: {data_changed}")

if attenuation_low < 0.1 and attenuation_high < 0.1 and preservation_pass > 0.5:
    print("\n✓ FILTER WORKING CORRECTLY!")
    print("  - Low frequencies attenuated by >10 dB")
    print("  - High frequencies attenuated by >10 dB")
    print("  - Passband preserved (>-3 dB)")
else:
    print("\n❌ FILTER MAY NOT BE WORKING AS EXPECTED")
    if attenuation_low >= 0.1:
        print(f"  - Low frequencies not sufficiently attenuated ({10*np.log10(attenuation_low):.1f} dB)")
    if attenuation_high >= 0.1:
        print(f"  - High frequencies not sufficiently attenuated ({10*np.log10(attenuation_high):.1f} dB)")
    if preservation_pass <= 0.5:
        print(f"  - Passband excessively attenuated ({10*np.log10(preservation_pass):.1f} dB)")
