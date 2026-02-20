"""
Test script for adaptive noise filtering

This script demonstrates how to use the adaptive_noise_filter method
to remove sensor cross-coupling interference that "walks" through seismic data.

The method extracts a noise template from a reference trace/window and uses
LMS adaptive filtering to subtract correlated noise from all traces.
"""

import sys
sys.path.insert(0, r'C:\Users\nieboer\gitclones\DeltaSEIS')

from deltaseis.readers.parser_seismic_segy import read_segy
import matplotlib.pyplot as plt
import numpy as np

# Configuration
CONFIG = {
    'data_path': r'C:\Users\nieboer\OneDrive - Stichting Deltares\data\xstar_heave\stacked',
    'segy_file': 'XStar_20241107_Heave-corrected_run-stacking_stacked_100_to_end.sgy',
    
    # Adaptive filter parameters
    'reference_trace': 1000,  # Trace with strong interference (adjust based on your data)
    'reference_start_ms': 50.0,  # Start of noise window (adjust to your interference location)
    'reference_end_ms': 55.0,   # End of noise window
    'filter_length': 100,       # Number of filter coefficients (~2 ms at 50 kHz)
    'max_lag_ms': 2.0,          # Maximum time shift to search for
    'mu': 0.01,                 # LMS step size (smaller = more stable, larger = faster)
    
    # Optional: apply only to specific trace range
    'apply_to_traces': None,  # None for all traces, or (start, end) tuple
    # 'apply_to_traces': (0, 5000),  # Example: only process first 5000 traces
    
    # Output
    'output_file': 'XStar_adaptive_filtered.sgy',
}

def main():
    print("="*70)
    print("Adaptive Noise Filtering Test")
    print("="*70)
    
    # Load data
    print("\n1. Loading seismic data...")
    segy_path = CONFIG['data_path'] + '\\' + CONFIG['segy_file']
    seismic = read_segy(segy_path)
    
    print(f"   Data shape: {seismic.data.shape}")
    print(f"   Sampling rate: {seismic.fs} Hz ({1e6/seismic.fs:.1f} μs interval)")
    print(f"   Number of traces: {seismic.data.shape[1]}")
    print(f"   Record length: {seismic.data.shape[0]/seismic.fs*1000:.2f} ms")
    
    # Show reference trace/window
    print("\n2. Visualizing reference noise template...")
    ref_trace = CONFIG['reference_trace']
    ref_start = CONFIG['reference_start_ms']
    ref_end = CONFIG['reference_end_ms']
    
    start_samp = int(ref_start * seismic.fs / 1000.0)
    end_samp = int(ref_end * seismic.fs / 1000.0)
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot reference trace
    ax = axes[0]
    time_vec = np.arange(seismic.data.shape[0]) / seismic.fs * 1000
    ax.plot(seismic.data[:, ref_trace], time_vec, 'k-', linewidth=0.5)
    ax.axhline(ref_start, color='r', linestyle='--', linewidth=2, label='Noise window')
    ax.axhline(ref_end, color='r', linestyle='--', linewidth=2)
    ax.fill_betweenx([ref_start, ref_end], ax.get_xlim()[0], ax.get_xlim()[1], 
                      color='red', alpha=0.2)
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Amplitude')
    ax.set_title(f'Reference Trace {ref_trace}')
    ax.invert_yaxis()
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot noise template
    ax = axes[1]
    noise_template = seismic.data[start_samp:end_samp, ref_trace]
    template_time = np.arange(len(noise_template)) / seismic.fs * 1000 + ref_start
    ax.plot(template_time, noise_template, 'r-', linewidth=1.5)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Extracted Noise Template')
    ax.grid(True, alpha=0.3)
    
    # Plot template spectrum
    ax = axes[2]
    fft_template = np.fft.rfft(noise_template)
    freqs = np.fft.rfftfreq(len(noise_template), 1.0/seismic.fs)
    ax.semilogy(freqs/1000, np.abs(fft_template)**2, 'r-', linewidth=1.5)
    ax.set_xlabel('Frequency (kHz)')
    ax.set_ylabel('Power Spectrum')
    ax.set_title('Noise Template Spectrum')
    ax.set_xlim([0, 20])
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('adaptive_filter_reference.png', dpi=150, bbox_inches='tight')
    print("   Saved: adaptive_filter_reference.png")
    plt.show()
    
    # Save original data for comparison
    print("\n3. Saving original data section for comparison...")
    data_original = seismic.data.copy()
    
    # Select a display window that includes some interference
    display_traces = slice(max(0, ref_trace - 50), min(seismic.data.shape[1], ref_trace + 50))
    display_time_ms = (0, 65)  # Full record length
    
    # Apply adaptive filtering
    print("\n4. Applying adaptive noise filter...")
    seismic.adaptive_noise_filter(
        reference_trace=CONFIG['reference_trace'],
        reference_start_ms=CONFIG['reference_start_ms'],
        reference_end_ms=CONFIG['reference_end_ms'],
        filter_length=CONFIG['filter_length'],
        max_lag_ms=CONFIG['max_lag_ms'],
        mu=CONFIG['mu'],
        apply_to_traces=CONFIG['apply_to_traces'],
        inplace=True
    )
    
    # Calculate noise reduction
    print("\n5. Calculating noise reduction metrics...")
    
    # Compute difference (removed noise)
    noise_removed = data_original - seismic.data
    
    # Calculate RMS values
    rms_original = np.sqrt(np.mean(data_original**2))
    rms_filtered = np.sqrt(np.mean(seismic.data**2))
    rms_noise = np.sqrt(np.mean(noise_removed**2))
    
    noise_reduction_db = 20 * np.log10(rms_original / (rms_filtered + 1e-10))
    noise_power_removed = (rms_noise / rms_original) * 100
    
    print(f"\n   RMS original: {rms_original:.6f}")
    print(f"   RMS filtered: {rms_filtered:.6f}")
    print(f"   RMS noise removed: {rms_noise:.6f}")
    print(f"   Noise reduction: {noise_reduction_db:.1f} dB")
    print(f"   Noise power removed: {noise_power_removed:.1f}%")
    
    # Visualization
    print("\n6. Creating comparison plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Calculate display parameters
    start_samp_display = int(display_time_ms[0] * seismic.fs / 1000.0)
    end_samp_display = int(display_time_ms[1] * seismic.fs / 1000.0)
    
    extent = [display_traces.start, display_traces.stop,
              display_time_ms[1], display_time_ms[0]]
    
    # Calculate clip value from original data
    clip_val = np.percentile(np.abs(data_original[start_samp_display:end_samp_display, display_traces]), 99)
    
    # Original data
    ax = axes[0, 0]
    im = ax.imshow(data_original[start_samp_display:end_samp_display, display_traces],
                   cmap='seismic', aspect='auto', extent=extent,
                   vmin=-clip_val, vmax=clip_val, interpolation='bilinear')
    ax.set_title('Original Data')
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Trace Number')
    # Mark reference trace and window
    ax.axvline(ref_trace, color='yellow', linestyle='--', linewidth=1, alpha=0.7)
    ax.axhline(ref_start, color='yellow', linestyle='--', linewidth=1, alpha=0.7)
    ax.axhline(ref_end, color='yellow', linestyle='--', linewidth=1, alpha=0.7)
    plt.colorbar(im, ax=ax, label='Amplitude')
    
    # Filtered data
    ax = axes[0, 1]
    im = ax.imshow(seismic.data[start_samp_display:end_samp_display, display_traces],
                   cmap='seismic', aspect='auto', extent=extent,
                   vmin=-clip_val, vmax=clip_val, interpolation='bilinear')
    ax.set_title('Adaptive Filtered Data')
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Trace Number')
    plt.colorbar(im, ax=ax, label='Amplitude')
    
    # Removed noise
    ax = axes[1, 0]
    noise_clip = np.percentile(np.abs(noise_removed[start_samp_display:end_samp_display, display_traces]), 99)
    im = ax.imshow(noise_removed[start_samp_display:end_samp_display, display_traces],
                   cmap='seismic', aspect='auto', extent=extent,
                   vmin=-noise_clip, vmax=noise_clip, interpolation='bilinear')
    ax.set_title(f'Removed Noise ({noise_power_removed:.1f}% of power)')
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Trace Number')
    plt.colorbar(im, ax=ax, label='Amplitude')
    
    # Single trace comparison
    ax = axes[1, 1]
    compare_trace = ref_trace  # Show the reference trace
    time_vec = np.arange(len(data_original[:, compare_trace])) / seismic.fs * 1000
    
    ax.plot(data_original[:, compare_trace], time_vec, 'r-', linewidth=1, alpha=0.7, label='Original')
    ax.plot(seismic.data[:, compare_trace], time_vec, 'b-', linewidth=1, label='Filtered')
    ax.plot(noise_removed[:, compare_trace] * 2, time_vec, 'g-', linewidth=0.5, alpha=0.7, label='Noise (×2)')
    
    ax.set_ylabel('Time (ms)')
    ax.set_xlabel('Amplitude')
    ax.set_title(f'Trace {compare_trace} Comparison')
    ax.invert_yaxis()
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('adaptive_filter_comparison.png', dpi=150, bbox_inches='tight')
    print("   Saved: adaptive_filter_comparison.png")
    plt.show()
    
    # Optional: Save filtered data
    print("\n7. Saving filtered data...")
    from deltaseis.export.export_seismic_segy import write_segy
    
    output_path = CONFIG['data_path'] + '\\' + CONFIG['output_file']
    write_segy(seismic, output_path, segy_path)
    print(f"   Saved: {CONFIG['output_file']}")
    
    print("\n" + "="*70)
    print("Adaptive filtering complete!")
    print("="*70)
    print("\nNext steps:")
    print("1. Examine 'adaptive_filter_comparison.png' to assess noise removal")
    print("2. If noise remains, try:")
    print("   - Adjusting reference_trace and reference_start/end_ms")
    print("   - Increasing filter_length for more complex noise patterns")
    print("   - Increasing max_lag_ms if noise has larger time delays")
    print("   - Adjusting mu (smaller for stability, larger for faster adaptation)")
    print("3. You can apply adaptive filtering multiple times with different")
    print("   reference templates if interference patterns vary")

if __name__ == '__main__':
    main()
