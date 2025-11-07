"""
Diagnostic script to analyze the signature wavelet extraction and frequency content
"""
import numpy as np
import matplotlib.pyplot as plt
import deltaseis
from pathlib import Path

# Load the data - use first file from the process_xstar.py segy_folder
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]

if len(segy_paths) == 0:
    raise FileNotFoundError(f"No SEG-Y files found in {segy_folder}")

segy_path = segy_paths[0]
print(f"Analyzing file: {segy_path.name}\n")

segy = deltaseis.Segy_edit(segy_path)
segy.xstar_split('full')

# Get trace array
trace_array = np.array(segy.trace_data)
print(f"Original data shape: {trace_array.shape}")

# Handle transpose if needed
if trace_array.shape[0] > trace_array.shape[1]:
    print(f"Transposing data from {trace_array.shape} to {trace_array.T.shape}")
    trace_array = trace_array.T

print(f"Final data shape (samples × traces): {trace_array.shape}")

# Initialize Seismic object
fs = 50_000  # 50 kHz sampling rate
seis = deltaseis.Seismic(trace_array, fs=fs, dx=0.4)

dt = 1/fs  # sampling interval
print(f"\nSeismic object dt: {dt*1e6:.2f} microseconds")
print(f"Total record length: {seis.data.shape[0] * dt * 1000:.2f} ms")
print(f"Number of traces: {seis.data.shape[1]}")

# Extract the wavelet manually to see what we're getting
trace_number = 6415  # User specified trace
start_time_ms = 30.5
end_time_ms = 31.7

# Convert to indices (note: trace_number is 1-indexed in user's mind, 0-indexed in code)
print(f"\n=== WAVELET EXTRACTION ===")
print(f"User requested trace: {trace_number}")
print(f"Python index will be: {trace_number - 1} (if 1-indexed) or {trace_number} (if 0-indexed)")

# Try both interpretations
for idx_offset, label in [(0, "0-indexed (trace as-is)"), (1, "1-indexed (trace - 1)")]:
    trace_idx = trace_number - idx_offset
    
    if trace_idx < 0 or trace_idx >= seis.data.shape[1]:
        print(f"\n{label}: INDEX OUT OF BOUNDS")
        continue
        
    print(f"\n{label}: Using trace index {trace_idx}")
    
    # Get time indices
    start_sample = int(start_time_ms * 1e-3 / dt)
    end_sample = int(end_time_ms * 1e-3 / dt)
    
    print(f"Time window: {start_time_ms} - {end_time_ms} ms")
    print(f"Sample indices: {start_sample} - {end_sample}")
    print(f"Wavelet length: {end_sample - start_sample} samples ({(end_sample - start_sample) * dt * 1000:.2f} ms)")
    
    # Extract wavelet
    signature = seis.data[start_sample:end_sample, trace_idx]
    print(f"Signature amplitude range: {np.min(signature):.2e} to {np.max(signature):.2e}")
    print(f"Signature RMS: {np.sqrt(np.mean(signature**2)):.2e}")
    
    # Compute frequency spectrum
    n_fft = seis.data.shape[0]  # Use same FFT length as full data
    sig_fft = np.fft.rfft(signature, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, dt)
    sig_power = np.abs(sig_fft)**2
    sig_amplitude = np.abs(sig_fft)
    
    # Find peak frequency
    peak_idx = np.argmax(sig_power)
    peak_freq = freqs[peak_idx]
    print(f"Peak frequency: {peak_freq:.0f} Hz")
    
    # Find frequencies above 50% of peak power
    half_power = np.max(sig_power) / 2
    strong_freqs = freqs[sig_power > half_power]
    if len(strong_freqs) > 0:
        print(f"Frequency range above half-power: {strong_freqs[0]:.0f} - {strong_freqs[-1]:.0f} Hz")
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Wavelet Analysis - {label}', fontsize=14, fontweight='bold')
    
    # Time domain - full trace
    time_full = np.arange(seis.data.shape[0]) * dt * 1000  # in ms
    axes[0, 0].plot(time_full, seis.data[:, trace_idx], 'b-', linewidth=0.5)
    axes[0, 0].axvline(start_time_ms, color='r', linestyle='--', label='Wavelet window')
    axes[0, 0].axvline(end_time_ms, color='r', linestyle='--')
    axes[0, 0].axhspan(np.min(signature), np.max(signature), 
                        xmin=start_time_ms/time_full[-1], xmax=end_time_ms/time_full[-1],
                        alpha=0.2, color='red')
    axes[0, 0].set_xlabel('Time (ms)')
    axes[0, 0].set_ylabel('Amplitude')
    axes[0, 0].set_title(f'Full Trace {trace_idx}')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Time domain - zoomed wavelet
    time_wavelet = np.arange(len(signature)) * dt * 1000 + start_time_ms
    axes[0, 1].plot(time_wavelet, signature, 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Time (ms)')
    axes[0, 1].set_ylabel('Amplitude')
    axes[0, 1].set_title(f'Extracted Wavelet ({len(signature)} samples, {(end_time_ms-start_time_ms):.1f} ms)')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Frequency spectrum - amplitude
    axes[1, 0].plot(freqs/1000, sig_amplitude, 'b-', linewidth=1)
    axes[1, 0].axvline(peak_freq/1000, color='r', linestyle='--', 
                       label=f'Peak: {peak_freq:.0f} Hz')
    axes[1, 0].set_xlabel('Frequency (kHz)')
    axes[1, 0].set_ylabel('Amplitude Spectrum')
    axes[1, 0].set_title('Frequency Content (Amplitude)')
    axes[1, 0].set_xlim([0, 25])
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Frequency spectrum - power (log scale)
    axes[1, 1].semilogy(freqs/1000, sig_power, 'g-', linewidth=1)
    axes[1, 1].axvline(peak_freq/1000, color='r', linestyle='--', 
                       label=f'Peak: {peak_freq:.0f} Hz')
    axes[1, 1].axhline(half_power, color='orange', linestyle='--', 
                       alpha=0.5, label='Half power')
    axes[1, 1].set_xlabel('Frequency (kHz)')
    axes[1, 1].set_ylabel('Power Spectrum (log scale)')
    axes[1, 1].set_title('Frequency Content (Power)')
    axes[1, 1].set_xlim([0, 25])
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_name = f'wavelet_analysis_trace{trace_idx}_{label.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")}.png'
    plt.savefig(output_name, dpi=150, bbox_inches='tight')
    print(f"\nSaved figure: {output_name}")
    
    plt.show()

print("\n=== ANALYSIS COMPLETE ===")
