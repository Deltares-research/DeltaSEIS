"""
Test and demonstrate time_frequency_denoise method
"""
import numpy as np
import matplotlib.pyplot as plt
import deltaseis
from pathlib import Path

# Load data
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]

if len(segy_paths) == 0:
    print("No data found!")
    exit()

segy_path = segy_paths[0]
print(f"Loading: {segy_path.name}\n")

s = deltaseis.Segy_edit(segy_path)
s.xstar_split('full')

# Work with subset for speed
trace_subset = np.array(s.trace_data).T[:, 5000:5300]
fs = 50_000
dt = 1/fs

print("="*70)
print("TIME-FREQUENCY DENOISING DEMONSTRATION")
print("="*70)

# Test 1: Original data
seis_orig = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_orig.time_power_gain(3)

# Test 2: STFT with soft thresholding (recommended)
seis_stft_soft = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_stft_soft.time_power_gain(3)
print("\n" + "-"*70)
print("TEST 1: STFT with Soft Thresholding (RECOMMENDED)")
print("-"*70)
seis_stft_soft.time_frequency_denoise(method='stft', threshold=0.1, 
                                      window_ms=10.0, threshold_type='soft')

# Test 3: STFT with hard thresholding
seis_stft_hard = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_stft_hard.time_power_gain(3)
print("\n" + "-"*70)
print("TEST 2: STFT with Hard Thresholding")
print("-"*70)
seis_stft_hard.time_frequency_denoise(method='stft', threshold=0.1, 
                                      window_ms=10.0, threshold_type='hard')

# Test 4: More aggressive denoising
seis_aggressive = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_aggressive.time_power_gain(3)
print("\n" + "-"*70)
print("TEST 3: Aggressive Denoising (threshold=0.2)")
print("-"*70)
seis_aggressive.time_frequency_denoise(method='stft', threshold=0.2, 
                                       window_ms=15.0, threshold_type='soft')

# Test 5: CWT denoising
seis_cwt = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis_cwt.time_power_gain(3)
print("\n" + "-"*70)
print("TEST 4: CWT Denoising (Multi-scale)")
print("-"*70)
try:
    seis_cwt.time_frequency_denoise(method='cwt', threshold=0.15, threshold_type='soft')
except ImportError:
    print("  PyWavelets not installed, skipping CWT test")
    seis_cwt = None

# Create comparison plots
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(3, 5, hspace=0.3, wspace=0.3)

time_ms = np.arange(trace_subset.shape[0]) * dt * 1000
extent = [0, trace_subset.shape[1], time_ms[-1], time_ms[0]]
vmax = np.percentile(np.abs(seis_orig.data), 98)

# Time domain comparisons
axes_time = [fig.add_subplot(gs[0, i]) for i in range(5)]
titles_time = ['Original', 'STFT Soft', 'STFT Hard', 'Aggressive', 'CWT']
datasets_time = [seis_orig, seis_stft_soft, seis_stft_hard, seis_aggressive, seis_cwt]

for ax, title, dataset in zip(axes_time, titles_time, datasets_time):
    if dataset is not None:
        im = ax.imshow(dataset.data, aspect='auto', cmap='seismic', 
                      vmin=-vmax, vmax=vmax, extent=extent)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('Trace')
        ax.set_ylabel('Time (ms)')
        ax.set_ylim([50, 20])
    else:
        ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=10)

# Frequency spectra
axes_freq = [fig.add_subplot(gs[1, i]) for i in range(5)]

freqs = np.fft.rfftfreq(trace_subset.shape[0], dt)

for ax, title, dataset in zip(axes_freq, titles_time, datasets_time):
    if dataset is not None:
        fft_data = np.fft.rfft(dataset.data, axis=0)
        power = np.mean(np.abs(fft_data)**2, axis=1)
        ax.semilogy(freqs/1000, power, linewidth=1.5)
        ax.set_xlabel('Frequency (kHz)')
        ax.set_ylabel('Power')
        ax.set_xlim([0, 20])
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{title} - Spectrum', fontsize=9)

# Single trace comparisons
axes_trace = [fig.add_subplot(gs[2, i]) for i in range(5)]
trace_idx = 150

for ax, title, dataset in zip(axes_trace, titles_time, datasets_time):
    if dataset is not None:
        ax.plot(dataset.data[:, trace_idx], time_ms, linewidth=1)
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Time (ms)')
        ax.set_ylim([50, 20])
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{title} - Trace {trace_idx}', fontsize=9)

fig.suptitle('Time-Frequency Denoising Comparison', fontsize=14, fontweight='bold', y=0.995)
plt.savefig('tf_denoise_comparison.png', dpi=150, bbox_inches='tight')
print("\n\nSaved: tf_denoise_comparison.png")
plt.show()

# Noise analysis
print("\n" + "="*70)
print("NOISE REDUCTION ANALYSIS")
print("="*70)

# Calculate noise estimate (difference from original)
def calculate_snr(signal, denoised):
    """Estimate SNR improvement."""
    noise_removed = signal - denoised
    signal_power = np.mean(signal**2)
    noise_power = np.mean(noise_removed**2)
    snr_db = 10 * np.log10(signal_power / (noise_power + 1e-10))
    return snr_db, noise_power

print("\nNoise removed (RMS of difference from original):")
for title, dataset in zip(titles_time[1:], datasets_time[1:]):
    if dataset is not None:
        snr_db, noise_power = calculate_snr(seis_orig.data, dataset.data)
        print(f"  {title:20s}: {np.sqrt(noise_power):.2e} (SNR: {snr_db:+.1f} dB)")

# Calculate coherency (trace-to-trace similarity)
def calculate_coherency(data):
    """Calculate average cross-correlation between adjacent traces."""
    n_samples, n_traces = data.shape
    coherency = []
    for i in range(n_traces - 1):
        corr = np.corrcoef(data[:, i], data[:, i+1])[0, 1]
        coherency.append(corr)
    return np.mean(coherency)

print("\nTrace-to-trace coherency (higher = more coherent signal):")
for title, dataset in zip(titles_time, datasets_time):
    if dataset is not None:
        coh = calculate_coherency(dataset.data)
        print(f"  {title:20s}: {coh:.4f}")

print("\n" + "="*70)
print("RECOMMENDATIONS")
print("="*70)
print("""
GENERAL GUIDELINES:
- Start with STFT + soft thresholding (fastest, most predictable)
- Typical threshold: 0.05-0.15 (lower = less aggressive)
- Window: 5-15 ms (shorter for better time resolution)
- Use CWT for signals with varying frequency content

PARAMETER SELECTION:
- If signal looks "washed out": decrease threshold (0.05-0.08)
- If noise still visible: increase threshold (0.15-0.25)
- If losing coherent events: increase window_ms (15-20 ms)
- If smearing events in time: decrease window_ms (5-8 ms)

WORKFLOW INTEGRATION:
1. Apply gain (time_power_gain)
2. Apply deconvolution (signature_deconvolution)
3. Apply TF denoising (removes amplified noise from decon)
4. Apply bandpass filter (removes remaining frequency-specific noise)
5. Apply top mute (removes water column)

QUALITY CONTROL:
- Check that coherent events are preserved
- Verify noise in mute zone is reduced
- Compare frequency spectra before/after
- Test on subset before processing full dataset
""")

print("\n" + "="*70)
print("USAGE EXAMPLES FOR process_xstar.py")
print("="*70)
print("""
# Example 1: Basic denoising (recommended starting point)
seis.time_frequency_denoise(threshold=0.1, window_ms=10.0)

# Example 2: Conservative denoising (preserve more signal)
seis.time_frequency_denoise(threshold=0.08, window_ms=12.0, threshold_type='soft')

# Example 3: Aggressive denoising (remove more noise)
seis.time_frequency_denoise(threshold=0.18, window_ms=15.0)

# Example 4: For data with varying noise levels
seis.time_frequency_denoise(threshold=0.12, global_threshold=False)

# Example 5: Using wavelets (slower but better for complex signals)
seis.time_frequency_denoise(method='cwt', threshold=0.15)

# Full workflow example:
seis.time_power_gain(3)
seis.signature_deconvolution(6414, 30.5, 31.7, epsilon=0.05)
seis.time_frequency_denoise(threshold=0.1, window_ms=10.0)  # ← Add this
seis.bandpass_filter(lowcut=2700, highcut=8000)
seis.top_mute(s.seabed_pick_savgol, shift_ms=-1.0, taper_ms=2.0)
""")
