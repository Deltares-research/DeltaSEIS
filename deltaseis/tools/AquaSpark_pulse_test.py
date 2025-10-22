#!/usr/bin/env python3
"""
AquaSpark Pulse Analyzer - Standard Version

Analyzes marine sparker pulse data from oscilloscope CSV recordings.
Processes multiple files, applies calibrations, and generates publication-quality plots.

Features:
- Automatic CSV file detection and processing
- Hydrophone calibration and unit conversion
- Signal stacking and spectrum analysis
- Publication-ready plot generation
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Tuple, Optional
import warnings

data_path = Path(r'D:\Projects\DIS_Offshore\pulse testing\scope_data_20180319_stacked_sparker\1')

# Configuration
CONFIG = {
    'hydrophone_sensitivity': 4.4,     # V/bar
    'hydrophone_depth': 5.5,          # m below sparker tips  
    'sample_rate': 0.0025,            # ms
    'max_frequency': 4000,            # Hz
    'plot_length': 3,                 # ms around peak
    'fft_padding': 3,                 # FFT zero-padding factor
    'invert_polarity': True,          # Multiply by -1
    'test_name': '30 random shots on 200 J',
    'data_dir': data_path,           # Default to current directory
    'figure_size': (16, 12),          # inches
    'line_width': 0.5,
    'font_size': 16,
    'time_ylim_min': -1.3,            # Bar-m units y-axis minimum
    'time_ylim_max': 1.3,             # Bar-m units y-axis maximum
    'spectra_ylim_min': -25,          # dB y-axis minimum
    'spectra_ylim_max': 5,            # dB y-axis maximum
}


def load_oscilloscope_data(filepath: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Load CSV data from oscilloscope with basic error handling."""
    try:
        # Skip header lines that start with # and read data
        # First, find where data starts
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        skip_rows = 0
        for i, line in enumerate(lines):
            if line.strip() and not line.strip().startswith('#'):
                skip_rows = i
                break
        
        # Read with pandas, skipping header lines
        df = pd.read_csv(filepath, skiprows=skip_rows, encoding='utf-8', on_bad_lines='skip')
        
        if df.shape[1] < 2:
            print(f"Warning: {filepath.name} has insufficient columns")
            return None, None
        
        time = df.iloc[:, 0].values.astype(float)
        amplitude = df.iloc[:, 1].values.astype(float)
        
        # Apply polarity inversion
        if CONFIG['invert_polarity']:
            amplitude *= -1
        
        return time, amplitude
        
    except Exception as e:
        print(f"Error loading {filepath.name}: {e}")
        return None, None


def process_trace(time: np.ndarray, amplitude: np.ndarray) -> Tuple[np.ndarray, ...]:
    """Process single trace: align, trim, and calibrate."""
    # Find peak and align time
    peak_idx = np.argmax(np.abs(amplitude))
    t_set = 0.003  # s
    time_aligned = time - time[peak_idx] + t_set
    
    # Trim around peak
    margin = 1550
    start_idx = max(0, peak_idx - margin)
    end_idx = min(len(amplitude), start_idx + 8000)
    
    time_trimmed = time_aligned[start_idx:end_idx]
    amp_trimmed = amplitude[start_idx:end_idx]
    
    # Convert to bar-m units
    sensitivity = CONFIG['hydrophone_sensitivity']
    depth = CONFIG['hydrophone_depth']
    bar_m = amp_trimmed * depth / sensitivity
    
    # Normalize
    normalized = amp_trimmed / np.max(np.abs(amp_trimmed))
    
    return time_trimmed, amp_trimmed, bar_m, normalized


def compute_spectrum(signal: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute power spectrum of signal."""
    # Calculate sampling frequency
    fs = 1000.0 / CONFIG['sample_rate']  # Hz
    
    # FFT with zero padding
    signal_length = len(signal)
    n_fft = 2 ** (int(np.log2(signal_length - 0.5)) + 1)
    padded_length = CONFIG['fft_padding'] * n_fft
    
    # Remove DC and compute FFT
    signal_dc = signal - np.mean(signal)
    fft_result = np.fft.fft(signal_dc, n=padded_length)
    
    # One-sided spectrum
    n_half = len(fft_result) // 2
    fft_half = fft_result[:n_half]
    
    # Convert to dB
    magnitude = np.abs(fft_half)
    power_db = 10 * np.log10(magnitude / np.max(magnitude))
    
    # Frequency axis
    frequencies = np.linspace(0, fs/2, n_half)
    
    return frequencies, power_db


def create_plots() -> Tuple[plt.Figure, tuple]:
    """Create figure with subplot layout."""
    # Configure matplotlib
    plt.rcParams.update({
        'font.size': CONFIG['font_size'],
        'figure.max_open_warning': 0,
    })
    
    # Create figure
    fig = plt.figure(figsize=CONFIG['figure_size'], dpi=100)
    fig.subplots_adjust(hspace=0.25, wspace=0.3, top=0.9, bottom=0.1, left=0.1, right=0.95)
    
    # Title
    title = f"Sparker Pulse Analysis: {CONFIG['test_name']}"
    if CONFIG['invert_polarity']:
        title += " (Inverted)"
    fig.suptitle(title, fontsize=CONFIG['font_size'] + 6, y=0.95)
    
    # Create subplots
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)  # Raw traces
    ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)  # Normalized traces  
    ax3 = plt.subplot2grid((2, 4), (1, 0), colspan=4)  # Spectra
    
    return fig, (ax1, ax2, ax3)


def setup_axes(axes: tuple) -> None:
    """Configure axis properties."""
    ax1, ax2, ax3 = axes
    
    # Raw traces (left)
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Output (Bar-m)')
    ax1.set_ylim(CONFIG['time_ylim_min'], CONFIG['time_ylim_max'])
    ax1.grid(visible=True, color='grey', linestyle='--', alpha=0.7)
    
    # Normalized traces (right)
    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Normalized (-)')
    ax2.set_ylim(-1.1, 1.1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.grid(visible=True, color='grey', linestyle='--', alpha=0.7)
    
    # Power spectra (bottom)
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylabel('Power (dB)')
    ax3.set_ylim(CONFIG['spectra_ylim_min'], CONFIG['spectra_ylim_max'])
    ax3.set_xlim(0, CONFIG['max_frequency'])
    ax3.grid(visible=True, color='grey', linestyle='--', alpha=0.7)


def main():
    """Main analysis pipeline."""
    # Find CSV files
    csv_files = list(CONFIG['data_dir'].glob('*.csv'))
    if not csv_files:
        print(f"No CSV files found in {CONFIG['data_dir']}")
        return
    
    print(f"Processing {len(csv_files)} files...")
    
    # Initialize storage
    traces_raw = []
    traces_bar = []
    traces_norm = []
    
    # Create plots
    fig, axes = create_plots()
    ax1, ax2, ax3 = axes
    
    # Generate colors
    colors = plt.cm.autumn(np.linspace(0, 1, len(csv_files)))
    
    # Process each file
    for i, filepath in enumerate(csv_files):
        time, amplitude = load_oscilloscope_data(filepath)
        if time is None:
            continue
            
        time_proc, amp_proc, bar_proc, norm_proc = process_trace(time, amplitude)
        
        # Store for stacking
        traces_raw.append(amp_proc)
        traces_bar.append(bar_proc)
        traces_norm.append(norm_proc)
        
        # Convert time to ms and set plot window
        time_ms = time_proc * 1000
        peak_idx = np.argmax(np.abs(amp_proc))
        t_min = time_ms[peak_idx] - 0.5
        t_max = t_min + CONFIG['plot_length']
        
        # Plot individual traces
        color = colors[i]
        label = filepath.stem
        
        ax1.plot(time_ms, bar_proc, color=color, 
                linewidth=CONFIG['line_width'], label=label, alpha=0.7)
        ax2.plot(time_ms, norm_proc, color=color, 
                linewidth=CONFIG['line_width'], label=label, alpha=0.7)
        
        # Compute and plot spectrum
        frequencies, power_db = compute_spectrum(amp_proc)
        ax3.plot(frequencies, power_db, color=color, 
                linewidth=CONFIG['line_width'], label=label, alpha=0.7)
        
        # Set time limits (use last trace for reference)
        ax1.set_xlim(t_min, t_max)
        ax2.set_xlim(t_min, t_max)
        
        print(f"Processed {filepath.name}: {1/(time[100]-time[99]):.0f} Hz, {len(time)} samples")
    
    # Compute and plot stacks
    if traces_bar:
        stack_bar = np.mean(traces_bar, axis=0)
        stack_norm = np.mean(traces_norm, axis=0) / np.max(np.abs(np.mean(traces_norm, axis=0)))
        
        # Plot stacks with emphasis
        ax1.plot(time_ms, stack_bar, 'blue', linewidth=CONFIG['line_width']*4, 
                label='Stack', zorder=10)
        ax2.plot(time_ms, stack_norm, 'blue', linewidth=CONFIG['line_width']*4, 
                label='Stack', zorder=10)
        
        # Stack spectrum
        stack_raw = np.mean(traces_raw, axis=0)
        freq_stack, power_stack = compute_spectrum(stack_raw)
        ax3.plot(freq_stack, power_stack, 'blue', linewidth=CONFIG['line_width']*4, 
                label='Stack', zorder=10)
        
        print(f"Stacked {len(traces_bar)} traces")
    
    # Finalize plot
    setup_axes(axes)
    
    # Save and show
    output_path = CONFIG['data_dir'] / f"stacked_{CONFIG['test_name'].replace(' ', '_')}.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    
    plt.tight_layout()
    plt.show()
    print("Analysis complete!")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Suppress font warnings
        main()