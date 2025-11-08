"""
Test the top_mute functionality with different parameters
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

# Get seabed pick
s.get_seabed_pick(10, 100, 9, 3, truncate=10)
s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)

print(f"\nSeabed pick statistics:")
print(f"  Min: {np.nanmin(s.seabed_pick):.2f} ms")
print(f"  Max: {np.nanmax(s.seabed_pick):.2f} ms")
print(f"  Mean: {np.nanmean(s.seabed_pick):.2f} ms")
print(f"  NaN count: {np.sum(np.isnan(s.seabed_pick))}")

# Work with subset for speed
trace_subset = np.array(s.trace_data).T[:, 5000:5500]
seabed_subset = s.seabed_pick[5000:5500]

fs = 50_000
dt = 1/fs

# Create comparison figure
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Top Mute Demonstration', fontsize=14, fontweight='bold')

time_ms = np.arange(trace_subset.shape[0]) * dt * 1000
extent = [0, trace_subset.shape[1], time_ms[-1], time_ms[0]]

# 1. Original data
seis_orig = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_orig.time_power_gain(3)

ax = axes[0, 0]
vmax = np.percentile(np.abs(seis_orig.data), 98)
im = ax.imshow(seis_orig.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset, 'g-', linewidth=2, label='Seabed pick')
ax.set_title('Original Data')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])
ax.legend()

# 2. Mute at seabed (no shift, no taper)
seis1 = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis1.time_power_gain(3)
seis1.top_mute(seabed_subset, shift_ms=0.0, taper_ms=0.0)

ax = axes[0, 1]
im = ax.imshow(seis1.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset, 'g-', linewidth=2, label='Mute boundary')
ax.set_title('Mute at Seabed (shift=0, no taper)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])
ax.legend()

# 3. Mute 1 ms above seabed (negative shift)
seis2 = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis2.time_power_gain(3)
seis2.top_mute(seabed_subset, shift_ms=-1.0, taper_ms=0.0)

ax = axes[0, 2]
im = ax.imshow(seis2.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset, 'g--', linewidth=1, alpha=0.5, label='Seabed pick')
ax.plot(np.arange(len(seabed_subset)), seabed_subset - 1.0, 'r-', linewidth=2, label='Mute boundary (-1 ms)')
ax.set_title('Mute 1 ms ABOVE Seabed (shift=-1 ms)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])
ax.legend()

# 4. Mute with 2 ms taper
seis3 = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis3.time_power_gain(3)
seis3.top_mute(seabed_subset, shift_ms=-1.0, taper_ms=2.0)

ax = axes[1, 0]
im = ax.imshow(seis3.data, aspect='auto', cmap='seismic', vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset - 1.0, 'r-', linewidth=2, label='Mute boundary')
# Show taper zone
taper_zone = ax.fill_between(np.arange(len(seabed_subset)), 
                              seabed_subset - 1.0, 
                              seabed_subset - 1.0 + 2.0,
                              alpha=0.3, color='yellow', label='Taper zone (2 ms)')
ax.set_title('With 2 ms Cosine Taper (shift=-1 ms)')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])
ax.legend()

# 5. Single trace comparison
trace_idx = 250
ax = axes[1, 1]
ax.plot(seis_orig.data[:, trace_idx], time_ms, 'k-', linewidth=1, alpha=0.5, label='Original')
ax.plot(seis1.data[:, trace_idx], time_ms, 'b-', linewidth=1.5, label='No taper')
ax.plot(seis3.data[:, trace_idx], time_ms, 'r-', linewidth=1.5, label='With taper')
ax.axhline(seabed_subset[trace_idx] - 1.0, color='g', linestyle='--', label='Mute boundary')
ax.axhline(seabed_subset[trace_idx] + 1.0, color='orange', linestyle='--', alpha=0.5, label='Taper end')
ax.set_xlabel('Amplitude')
ax.set_ylabel('Time (ms)')
ax.set_title(f'Single Trace Comparison (trace {trace_idx})')
ax.set_ylim([50, 20])
ax.legend()
ax.grid(True, alpha=0.3)

# 6. Amplitude comparison in taper zone
ax = axes[1, 2]
# Extract a window around the taper zone
taper_center_sample = int((seabed_subset[trace_idx] - 1.0 + 1.0) * 1e-3 / dt)
window = slice(max(0, taper_center_sample - 150), min(len(time_ms), taper_center_sample + 150))
ax.plot(np.abs(seis_orig.data[window, trace_idx]), time_ms[window], 'k-', linewidth=1, alpha=0.5, label='Original')
ax.plot(np.abs(seis1.data[window, trace_idx]), time_ms[window], 'b-', linewidth=1.5, label='Hard cutoff')
ax.plot(np.abs(seis3.data[window, trace_idx]), time_ms[window], 'r-', linewidth=1.5, label='Tapered')
ax.axhline(seabed_subset[trace_idx] - 1.0, color='g', linestyle='--', alpha=0.7)
ax.axhspan(seabed_subset[trace_idx] - 1.0, seabed_subset[trace_idx] + 1.0, 
           alpha=0.2, color='yellow')
ax.set_xlabel('Absolute Amplitude')
ax.set_ylabel('Time (ms)')
ax.set_title('Amplitude Envelope - Taper Effect')
ax.set_ylim([seabed_subset[trace_idx] - 5, seabed_subset[trace_idx] + 5])
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('top_mute_demonstration.png', dpi=150, bbox_inches='tight')
print("\nSaved: top_mute_demonstration.png")
plt.show()

print("\n" + "="*60)
print("USAGE EXAMPLES FOR process_xstar.py")
print("="*60)
print("\n# Example 1: Mute 1 ms above seabed, no taper")
print("seis.top_mute(s.seabed_pick, shift_ms=-1.0)")
print("\n# Example 2: Mute 2 ms above seabed with 3 ms taper")
print("seis.top_mute(s.seabed_pick, shift_ms=-2.0, taper_ms=3.0)")
print("\n# Example 3: Mute using filtered seabed pick")
print("seis.top_mute(s.seabed_pick_savgol, shift_ms=-1.5, taper_ms=2.0)")
print("\n# Example 4: Mute at constant time (e.g., 25 ms)")
print("seis.top_mute(25.0, taper_ms=2.0)")
