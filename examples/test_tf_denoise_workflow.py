"""
Practical example: TF denoise in a processing workflow
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
s.get_seabed_pick(10, 100, 9, 3, truncate=10)
s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)

# Work with subset
trace_subset = np.array(s.trace_data).T[:, 6000:6500]
seabed_subset = s.seabed_pick_savgol[6000:6500]
fs = 50_000
dt = 1/fs

print("="*70)
print("PROCESSING WORKFLOW WITH TF DENOISING")
print("="*70)

# Workflow WITHOUT TF denoising
print("\nWorkflow A: WITHOUT TF Denoising")
print("-"*70)
seis_without = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_without.time_power_gain(3)
seis_without.signature_deconvolution(415, 30.5, 31.7, method='wiener', 
                                    epsilon=0.05, prewhiten=True, prewhiten_percent=1.0)
seis_without.bandpass_filter(lowcut=2700, highcut=8000)
seis_without.top_mute(seabed_subset, shift_ms=-1.0, taper_ms=2.0)

# Workflow WITH TF denoising (RECOMMENDED)
print("\n" + "="*70)
print("Workflow B: WITH TF Denoising (After Decon)")
print("-"*70)
seis_with = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_with.time_power_gain(3)
seis_with.signature_deconvolution(415, 30.5, 31.7, method='wiener', 
                                  epsilon=0.05, prewhiten=True, prewhiten_percent=1.0)
# ↓↓↓ ADD TF DENOISING HERE ↓↓↓
seis_with.time_frequency_denoise(threshold=0.12, window_ms=10.0, threshold_type='soft')
# ↑↑↑ This removes noise amplified by deconvolution ↑↑↑
seis_with.bandpass_filter(lowcut=2700, highcut=8000)
seis_with.top_mute(seabed_subset, shift_ms=-1.0, taper_ms=2.0)

# Create comparison
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Impact of TF Denoising in Processing Workflow', fontsize=14, fontweight='bold')

time_ms = np.arange(trace_subset.shape[0]) * dt * 1000
extent = [0, trace_subset.shape[1], time_ms[-1], time_ms[0]]
vmax = np.percentile(np.abs(seis_without.data), 98)

# Row 1: Seismic sections
ax = axes[0, 0]
im = ax.imshow(seis_without.data, aspect='auto', cmap='seismic', 
               vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset, 'g-', linewidth=2, alpha=0.7)
ax.set_title('Without TF Denoising', fontsize=11, fontweight='bold')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

ax = axes[0, 1]
im = ax.imshow(seis_with.data, aspect='auto', cmap='seismic', 
               vmin=-vmax, vmax=vmax, extent=extent)
ax.plot(np.arange(len(seabed_subset)), seabed_subset, 'g-', linewidth=2, alpha=0.7)
ax.set_title('With TF Denoising', fontsize=11, fontweight='bold')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

# Difference (noise removed)
ax = axes[0, 2]
diff = seis_without.data - seis_with.data
vmax_diff = np.percentile(np.abs(diff), 98)
im = ax.imshow(diff, aspect='auto', cmap='seismic', 
               vmin=-vmax_diff, vmax=vmax_diff, extent=extent)
ax.set_title('Noise Removed by TF Denoising', fontsize=11, fontweight='bold')
ax.set_xlabel('Trace')
ax.set_ylabel('Time (ms)')
ax.set_ylim([50, 20])

# Row 2: Frequency spectra
freqs = np.fft.rfftfreq(trace_subset.shape[0], dt)
fft_without = np.fft.rfft(seis_without.data, axis=0)
fft_with = np.fft.rfft(seis_with.data, axis=0)
power_without = np.mean(np.abs(fft_without)**2, axis=1)
power_with = np.mean(np.abs(fft_with)**2, axis=1)

ax = axes[1, 0]
ax.semilogy(freqs/1000, power_without, 'b-', linewidth=2, label='Without TF denoise')
ax.semilogy(freqs/1000, power_with, 'r-', linewidth=2, label='With TF denoise')
ax.axvline(2.7, color='green', linestyle='--', alpha=0.5, label='Bandpass range')
ax.axvline(8.0, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('Frequency (kHz)')
ax.set_ylabel('Average Power')
ax.set_xlim([0, 20])
ax.set_title('Frequency Content Comparison', fontsize=11)
ax.legend()
ax.grid(True, alpha=0.3)

# Single trace comparison
ax = axes[1, 1]
trace_idx = 250
ax.plot(seis_without.data[:, trace_idx], time_ms, 'b-', 
        linewidth=1.5, alpha=0.7, label='Without')
ax.plot(seis_with.data[:, trace_idx], time_ms, 'r-', 
        linewidth=1.5, label='With TF denoise')
ax.axhline(seabed_subset[trace_idx], color='g', linestyle='--', alpha=0.5)
ax.set_xlabel('Amplitude')
ax.set_ylabel('Time (ms)')
ax.set_title(f'Single Trace (#{trace_idx})', fontsize=11)
ax.set_ylim([50, 20])
ax.legend()
ax.grid(True, alpha=0.3)

# Zoomed view around seabed
ax = axes[1, 2]
zoom_center = seabed_subset[trace_idx]
zoom_window = 10  # ms
zoom_idx_start = int((zoom_center - zoom_window) * 1e-3 / dt)
zoom_idx_end = int((zoom_center + zoom_window) * 1e-3 / dt)
zoom_idx_start = max(0, zoom_idx_start)
zoom_idx_end = min(len(time_ms), zoom_idx_end)

ax.plot(seis_without.data[zoom_idx_start:zoom_idx_end, trace_idx], 
        time_ms[zoom_idx_start:zoom_idx_end], 'b-', 
        linewidth=1.5, alpha=0.7, label='Without')
ax.plot(seis_with.data[zoom_idx_start:zoom_idx_end, trace_idx], 
        time_ms[zoom_idx_start:zoom_idx_end], 'r-', 
        linewidth=1.5, label='With TF denoise')
ax.axhline(seabed_subset[trace_idx], color='g', linestyle='--', 
           linewidth=2, alpha=0.7, label='Seabed')
ax.set_xlabel('Amplitude')
ax.set_ylabel('Time (ms)')
ax.set_title(f'Zoomed View (±{zoom_window} ms)', fontsize=11)
ax.invert_yaxis()
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('tf_denoise_workflow_comparison.png', dpi=150, bbox_inches='tight')
print("\n\nSaved: tf_denoise_workflow_comparison.png")
plt.show()

# Statistics
print("\n" + "="*70)
print("QUANTITATIVE COMPARISON")
print("="*70)

noise_removed = seis_without.data - seis_with.data
noise_rms = np.sqrt(np.mean(noise_removed**2))
signal_rms_without = np.sqrt(np.mean(seis_without.data**2))
signal_rms_with = np.sqrt(np.mean(seis_with.data**2))

print(f"\nRMS Amplitude:")
print(f"  Without TF denoising: {signal_rms_without:.4f}")
print(f"  With TF denoising:    {signal_rms_with:.4f}")
print(f"  Noise removed (RMS):  {noise_rms:.4f}")
print(f"  Noise fraction:       {(noise_rms/signal_rms_without)*100:.1f}%")

# Check data in mute zone
mute_zone_samples = []
for i, t in enumerate(seabed_subset):
    mute_sample = int((t - 1.0) * 1e-3 / dt)  # 1 ms above seabed
    if mute_sample > 0:
        mute_zone_samples.append(seis_without.data[:mute_sample, i])
        mute_zone_samples.append(seis_with.data[:mute_sample, i])

if len(mute_zone_samples) > 0:
    mute_rms_without = np.sqrt(np.mean(np.concatenate([seis_without.data[:int(np.nanmin(seabed_subset)*1e-3/dt), :]])**2))
    mute_rms_with = np.sqrt(np.mean(np.concatenate([seis_with.data[:int(np.nanmin(seabed_subset)*1e-3/dt), :]])**2))
    
    print(f"\nNoise in mute zone (above seabed):")
    print(f"  Without TF denoising: {mute_rms_without:.4f}")
    print(f"  With TF denoising:    {mute_rms_with:.4f}")
    print(f"  Reduction:            {(1 - mute_rms_with/mute_rms_without)*100:.1f}%")

print("\n" + "="*70)
print("RECOMMENDATION")
print("="*70)
print("""
✓ TF denoising is MOST EFFECTIVE when applied AFTER deconvolution:
  
  Reason: Deconvolution amplifies both signal AND noise.
          TF denoising removes the amplified random noise while
          preserving the coherent signal.

Processing order:
  1. time_power_gain()           ← Apply gain
  2. signature_deconvolution()   ← Improves resolution (amplifies noise!)
  3. time_frequency_denoise()    ← Remove amplified noise ✓✓✓
  4. bandpass_filter()           ← Remove out-of-band noise
  5. top_mute()                  ← Remove water column

Typical parameters for XStar data:
  threshold=0.10-0.15  (start conservative, increase if needed)
  window_ms=10.0       (good balance for 2-8 kHz data)
  threshold_type='soft' (smoother results)
""")
