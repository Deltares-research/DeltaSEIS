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

data_path = Path(r'c:\Users\nieboer\OneDrive - Stichting Deltares\pulse testing\dis2025\Diepte38cm')

# Configuration
CONFIG = {
    'hydrophone_sensitivity': 4.4,    # V/bar
    'hydrophone_depth': 5.5,          # m below sparker tips  
    'sample_rate': 6250,              # kHz
    'max_plot_frequency': 8000,       # Hz
    'plot_length': 4.0,                 # ms around peak
    'fft_padding': 3,                 # FFT zero-padding factor
    'invert_polarity': True,          # Multiply by -1
    'test_name': '30 random shots on 200 J',
    'data_dir': data_path,           # Default to current directory
    'figure_size': (12, 8),          # inches
    'line_width': 0.5,
    'font_size': 16,
    'time_ylim_min': -1.3,            # Bar-m units y-axis minimum
    'time_ylim_max': 1.3,             # Bar-m units y-axis maximum
    'spectra_ylim_min': -15,          # dB y-axis minimum
    'spectra_ylim_max': 5,            # dB y-axis maximum,
    'alignment_time': 3.0,            # ms, time to align peak to
    'start_time': 2.0,               # ms, start window (should be < alignment_time)
    'end_time': 4.0,                  # ms, end window (should be > alignment_time, set to high number to include full trace)
    'number_of_best_shots': 100,        # Number of best shots to select by amplitude
    'subfolder_name_filter': 'hz'      # Set to None to select all subfolders, or a string to filter
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
    
    # Title will be set to figure file name later in main()
    
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
    ax1.set_ylim(-6, 12)
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
    """Main analysis pipeline: iterate over subfolders and save figures in figures/ subfolder."""
    data_dir = CONFIG['data_dir']
    figures_dir = data_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    subfolder_name_filter = CONFIG.get('subfolder_name_filter', None)
    subfolders = [f for f in data_dir.iterdir() if f.is_dir() and (subfolder_name_filter is None or subfolder_name_filter.lower() in f.name.lower())]
    if not subfolders:
        print(f"No subfolders found in {data_dir}")
        return


    for subfolder in subfolders:
        csv_files = list(subfolder.glob('*.csv'))
        if not csv_files:
            print(f"No CSV files found in {subfolder}")
            continue

        print(f"Processing {len(csv_files)} files in {subfolder.name}...")

        # Load and process all traces
        traces_raw = []
        traces_bar = []
        traces_norm = []
        amplitudes = []
        time_ms_list = []
        sample_counts = []
        for filepath in csv_files:
            time, amplitude = load_oscilloscope_data(filepath)
            if time is None:
                continue
            time_proc, amp_proc, bar_proc, norm_proc = process_trace(time, amplitude)
            traces_raw.append(amp_proc)
            traces_bar.append(bar_proc)
            traces_norm.append(norm_proc)
            amplitudes.append(np.max(np.abs(amp_proc)))
            time_ms_list.append(time_proc * 1000)
            sample_counts.append(len(amp_proc))

        # Select best shots by amplitude
        n_best = min(CONFIG['number_of_best_shots'], len(traces_bar))
        if CONFIG['number_of_best_shots'] > len(traces_bar):
            print(f"Warning: Requested {CONFIG['number_of_best_shots']} best shots, but only {len(traces_bar)} files available. Using all files.")
        best_indices = np.argsort(amplitudes)[-n_best:][::-1]  # descending order
        traces_raw = [traces_raw[i] for i in best_indices]
        traces_bar = [traces_bar[i] for i in best_indices]
        traces_norm = [traces_norm[i] for i in best_indices]
        time_ms_list = [time_ms_list[i] for i in best_indices]

        # Create plots
        fig, axes = create_plots()
        ax1, ax2, ax3 = axes

        # Generate colors
        colors = plt.cm.autumn(np.linspace(0, 1, len(traces_bar)))

        # Plot each best shot
        for i in range(len(traces_bar)):
            time_ms = time_ms_list[i]
            bar_proc = traces_bar[i]
            norm_proc = traces_norm[i]
            amp_proc = traces_raw[i]
            peak_idx = np.argmax(np.abs(amp_proc))
            t_min = time_ms[peak_idx] - 0.5
            t_max = t_min + CONFIG['plot_length']
            color = colors[i]
            ax1.plot(time_ms, bar_proc, color=color,
                     linewidth=CONFIG['line_width'], label=f"Shot {i+1}", alpha=0.7)
            ax2.plot(time_ms, norm_proc, color=color,
                     linewidth=CONFIG['line_width'], label=f"Shot {i+1}", alpha=0.7)
            frequencies, power_db = compute_spectrum(amp_proc)
            ax3.plot(frequencies, power_db, color=color, 
                    linewidth=CONFIG['line_width'], label=f"Shot {i+1}", alpha=0.7)
            # Adjust x-axis limits if CONFIG['end_time'] is less than t_max
            end_time_ms = CONFIG['end_time']

            if end_time_ms < t_max:
                ax1.set_xlim(t_min, end_time_ms)
                ax2.set_xlim(t_min, end_time_ms)
            else:
                ax1.set_xlim(t_min, t_max)
                ax2.set_xlim(t_min, t_max)
            
                ax1.set_xlim(t_min, t_max)
                ax2.set_xlim(t_min, t_max)

        # Set figure title to output file name
        output_path = figures_dir / f"{subfolder.name}_stacked.png"
        fig.suptitle(str(output_path.name), fontsize=CONFIG['font_size'] + 6, y=1.01)

        # Compute and plot stacks
        if traces_bar:
            # Find minimum length among all traces
            min_len = min(len(tr) for tr in traces_bar)
            # Trim all traces to min_len
            traces_bar_trimmed = [tr[:min_len] for tr in traces_bar]
            traces_norm_trimmed = [tr[:min_len] for tr in traces_norm]
            traces_raw_trimmed = [tr[:min_len] for tr in traces_raw]
            # Also trim time_ms for plotting
            time_ms_stack = time_ms[:min_len]

            stack_bar = np.mean(traces_bar_trimmed, axis=0)
            stack_norm = np.mean(traces_norm_trimmed, axis=0) / np.max(np.abs(np.mean(traces_norm_trimmed, axis=0)))

            # Plot stacks with emphasis
            ax1.plot(time_ms_stack, stack_bar, 'blue', linewidth=CONFIG['line_width']*4, 
                    label='Stack', zorder=10)
            ax2.plot(time_ms_stack, stack_norm, 'blue', linewidth=CONFIG['line_width']*4, 
                    label='Stack', zorder=10)

            # Stack spectrum
            stack_raw = np.mean(traces_raw_trimmed, axis=0)
            freq_stack, power_stack = compute_spectrum(stack_raw)
            ax3.plot(freq_stack, power_stack, 'blue', linewidth=CONFIG['line_width']*4, 
                    label='Stack', zorder=10)

            print(f"Stacked {len(traces_bar)} traces.")

            # Finalize plot
            setup_axes(axes)

            # Save figure in figures/ subfolder of data_path
            output_path = figures_dir / f"{subfolder.name}_stacked.png"
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {output_path}")

            plt.close(fig)

    print("Analysis complete!")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Suppress font warnings
        main()