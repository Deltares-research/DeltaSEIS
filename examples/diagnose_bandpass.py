"""
Diagnose bandpass filter issue
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
trace_subset = np.array(s.trace_data).T[:, 5000:5100]
fs = 50_000
dt = 1/fs

# Create Seismic object
seis = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis.time_power_gain(3)

print("="*60)
print("TESTING BANDPASS FILTER")
print("="*60)

# Check attributes before filter
print("\nBefore bandpass_filter:")
print(f"  hasattr(seis, 'data'): {hasattr(seis, 'data')}")
print(f"  hasattr(seis, 'data_bandpass'): {hasattr(seis, 'data_bandpass')}")
print(f"  seis.data.shape: {seis.data.shape}")

# Store copy of original
data_before = seis.data.copy()

# Apply bandpass filter
print("\nApplying bandpass_filter(2700, 8000)...")
seis.bandpass_filter(lowcut=2700, highcut=8000)

# Check attributes after filter
print("\nAfter bandpass_filter:")
print(f"  hasattr(seis, 'data'): {hasattr(seis, 'data')}")
print(f"  hasattr(seis, 'data_bandpass'): {hasattr(seis, 'data_bandpass')}")

if hasattr(seis, 'data_bandpass'):
    print(f"  seis.data_bandpass.shape: {seis.data_bandpass.shape}")
    print(f"  seis.data.shape: {seis.data.shape}")
    
    # Check if data was modified
    data_changed = not np.allclose(seis.data, data_before)
    print(f"\n  seis.data was modified: {data_changed}")
    
    data_bandpass_different = not np.allclose(seis.data_bandpass, data_before[:, 1:])
    print(f"  seis.data_bandpass is different from original: {data_bandpass_different}")

# Compute frequency spectra
fft_before = np.fft.rfft(data_before, axis=0)
fft_after = np.fft.rfft(seis.data, axis=0)
freqs = np.fft.rfftfreq(data_before.shape[0], dt)

power_before = np.mean(np.abs(fft_before)**2, axis=1)
power_after = np.mean(np.abs(fft_after)**2, axis=1)

# Plot comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Bandpass Filter Diagnostic', fontsize=14, fontweight='bold')

# Frequency spectrum
ax = axes[0, 0]
ax.semilogy(freqs/1000, power_before, 'b-', linewidth=1.5, label='Before filter')
ax.semilogy(freqs/1000, power_after, 'r-', linewidth=1.5, label='After filter (seis.data)')
if hasattr(seis, 'data_bandpass'):
    fft_bandpass = np.fft.rfft(seis.data_bandpass, axis=0)
    power_bandpass = np.mean(np.abs(fft_bandpass)**2, axis=1)
    ax.semilogy(freqs/1000, power_bandpass, 'g--', linewidth=2, label='data_bandpass attribute')
ax.axvline(2.7, color='orange', linestyle='--', alpha=0.7, label='Lowcut (2.7 kHz)')
ax.axvline(8.0, color='orange', linestyle='--', alpha=0.7, label='Highcut (8.0 kHz)')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Power')
ax.set_xlim([0, 15])
ax.set_title('Frequency Spectrum Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# Ratio plot
ax = axes[0, 1]
ratio = power_after / (power_before + 1e-10)
ax.plot(freqs/1000, ratio, 'r-', linewidth=1.5, label='After / Before')
ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
ax.axvline(2.7, color='orange', linestyle='--', alpha=0.7)
ax.axvline(8.0, color='orange', linestyle='--', alpha=0.7)
ax.axhspan(0, 1, xmin=0, xmax=2.7/15, alpha=0.2, color='red', label='Should be attenuated')
ax.axhspan(0, 1, xmin=8.0/15, xmax=1, alpha=0.2, color='red')
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Power Ratio (After/Before)')
ax.set_xlim([0, 15])
ax.set_ylim([0, 1.2])
ax.set_title('Filter Response (should be <1 outside passband)')
ax.legend()
ax.grid(True, alpha=0.3)

# Time domain - before
ax = axes[1, 0]
time_ms = np.arange(data_before.shape[0]) * dt * 1000
extent = [0, data_before.shape[1], time_ms[-1], time_ms[0]]
vmax = np.percentile(np.abs(data_before), 98)
im = ax.imshow(data_before, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('Before Filter')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

# Time domain - after
ax = axes[1, 1]
im = ax.imshow(seis.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.set_title('After Filter (seis.data)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

plt.tight_layout()
plt.savefig('bandpass_filter_diagnostic.png', dpi=150, bbox_inches='tight')
print("\nSaved: bandpass_filter_diagnostic.png")
plt.show()

print("\n" + "="*60)
print("DIAGNOSIS")
print("="*60)

if np.allclose(seis.data, data_before):
    print("\n❌ BUG CONFIRMED: seis.data was NOT modified by bandpass_filter!")
    print("   The filter creates 'data_bandpass' attribute instead of updating 'data'")
    print("\n   This is a bug in deltaseis/base_seismic.py line 358-361:")
    print("   Should be: self.data = data_bandpass")
    print("   Currently: self.data_bandpass = data_bandpass[:, 1:]")
else:
    print("\n✓ seis.data was modified by the filter")

if hasattr(seis, 'data_bandpass'):
    print("\n   Note: data_bandpass attribute exists but is not being used")
    print(f"   It has shape {seis.data_bandpass.shape} (first trace removed - another bug!)")
