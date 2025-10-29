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

data_path = Path(r'D:\Projects\DIS_Offshore\pulse testing\dis2025\400J')

# Configuration
CONFIG = {
    'hydrophone_sensitivity': 4.4,    # V/bar
    'hydrophone_depth': 5.5,          # m below sparker tips  
    'sample_rate': 6250,              # kHz
    'max_plot_frequency': 4000,       # Hz
    'plot_length': 3,                 # ms around peak
    'fft_padding': 3,                 # FFT zero-padding factor
    'invert_polarity': True,          # Multiply by -1
    'test_name': '30 random shots on 200 J',
    'data_dir': data_path,           # Default to current directory
    'figure_size': (12, 8),          # inches
    'line_width': 0.5,
    'font_size': 16,
    'time_ylim_min': -1.3,            # Bar-m units y-axis minimum
    'time_ylim_max': 1.3,             # Bar-m units y-axis maximum
    'spectra_ylim_min': -50,          # dB y-axis minimum
    'spectra_ylim_max': 5,            # dB y-axis maximum,
    'alignment_time': 3.0,            # ms, time to align peak to
    'start_time': 2.0,               # ms, start window (should be < alignment_time)
    'end_time': 6.0,                  # ms, end window (should be > alignment_time)
    'figure_title': "GSO 360 sparker pulse analysis (400 J)"
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
    t_set = CONFIG['alignment_time'] * 1e-3  # s
    tstart = CONFIG['start_time'] * 1e-3  # s
    tend = CONFIG['end_time'] * 1e-3  # s

    peak_idx = np.argmax(np.abs(amplitude))
    
    # change the time axis so that all peaks are aligned at t_set (from arbitrary user input alignment time)
    time_aligned = time - time[peak_idx] + t_set
    
    # Trim around peak
    start_idx = np.argmin(np.abs(time_aligned - tstart))
    end_idx = np.argmin(np.abs(tend - time_aligned))
    
       
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
    """Compute power spectrum of signal with proper windowing and scaling."""
    
    # Sampling frequency
    fs = CONFIG['sample_rate'] * 1e3  # Hz
    
    # Remove DC
    signal = signal - np.mean(signal)
    
    # Apply Hann window
    window = np.hanning(len(signal))
    signal_win = signal * window
    
    # FFT (one-sided)
    fft_vals = np.fft.rfft(signal_win)
    freqs = np.fft.rfftfreq(len(signal), 1/fs)
    
    # Amplitude scaling (compensate window energy)
    magnitude = np.abs(fft_vals) / (np.sum(window) / 2)
    
    # Convert to dBV (reference = 1 V)
    power_dbv = 20 * np.log10(magnitude + 1e-12)  # avoid log(0)
    
    return freqs, power_dbv


def create_plots() -> Tuple[plt.Figure, tuple]:
    """Create figure with subplot layout."""
    # Configure matplotlib
    plt.rcParams.update({
        'font.size': CONFIG['font_size'],
        'figure.max_open_warning': 0,
    })
    
    # Create figure
    fig = plt.figure(figsize=CONFIG['figure_size'], dpi=100)
    fig.subplots_adjust(hspace=0.25, wspace=0.3, top=0.95, bottom=0.07, left=0.08, right=0.9)
    
    # Title
    title = CONFIG['figure_title']
    if CONFIG['invert_polarity']:
        title += " (Inverted)"
    fig.suptitle(title, fontsize=CONFIG['font_size'] + 6, y=1.01)
    
    # Create subplots
    ax1 = plt.subplot2grid((5, 4), (0, 0), colspan=2, rowspan=3)  # Raw traces (60%)
    ax2 = plt.subplot2grid((5, 4), (0, 2), colspan=2, rowspan=3)  # Normalized traces (60%)
    ax3 = plt.subplot2grid((5, 4), (3, 0), colspan=4, rowspan=2)  # Spectra (40%)
    
    return fig, (ax1, ax2, ax3)


def setup_axes(axes: tuple) -> None:
    """Configure axis properties."""
    ax1, ax2, ax3 = axes
    
    # Raw traces (left)
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Output (Bar-m)')
    ax1.set_ylim(-4, 8)
    ax1.grid(visible=True, color='grey', linestyle='--', alpha=0.7)
    
    # Normalized traces (right)
    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Normalized (-)')
    ax2.set_ylim(-0.8, 1.2)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.grid(visible=True, color='grey', linestyle='--', alpha=0.7)
    
    # Power spectra (bottom)
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylabel('Power (dB)')
    ax3.set_ylim(CONFIG['spectra_ylim_min'], CONFIG['spectra_ylim_max'])
    ax3.set_xlim(0, CONFIG['max_plot_frequency'])
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
    plots_dir = Path(__file__).parent / "plots"
    plots_dir.mkdir(exist_ok=True)
    output_path = plots_dir / f"stacked_{CONFIG['test_name'].replace(' ', '_')}.png"
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    
    plt.tight_layout()
    plt.show()
    print("Analysis complete!")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Suppress font warnings
        main()