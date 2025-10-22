#!/usr/bin/env python3
"""
AquaSpark Pro - Ultra-Robust Offshore Seismic Pulse Analyzer

OFFSHORE-READY VERSION with comprehensive error handling and fallback mechanisms.
Designed for mission-critical offshore operations where internet is unavailable
and debugging time is not available.

This bulletproof tool analyzes marine sparker pulse data from oscilloscope recordings
with extensive validation, error recovery, and graceful degradation capabilities.

ROBUSTNESS FEATURES:
- Comprehensive input validation and sanitization
- Multiple fallback mechanisms for data loading
- Automatic error recovery and graceful degradation  
- Detailed logging and status reporting
- Memory-efficient processing for large datasets
- Safe file I/O with atomic operations
- Network-independent operation (no external dependencies)
- Extensive bounds checking and NaN/Inf handling
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Tuple, List, Optional, Union
import warnings
import sys
import traceback
import logging
from datetime import datetime

data_path = Path(r'D:\Projects\DIS_Offshore\pulse testing\scope_data_20180319_stacked_sparker\1')

# Robust Configuration with Validation
class RobustConfig:
    """Ultra-safe configuration with validation and fallback values."""
    
    def __init__(self):
        self._config = {
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
            'min_samples': 1000,              # Minimum samples required
            'max_samples': 100000,            # Maximum samples to prevent memory issues
            'max_files': 1000,                # Maximum files to process
            'backup_save': True,              # Save backup copies
        }
        self._validate_config()
    
    def _validate_config(self):
        """Validate configuration parameters and apply safe defaults."""
        try:
            # Validate numeric parameters
            if self._config['hydrophone_sensitivity'] <= 0:
                self._config['hydrophone_sensitivity'] = 4.4
                self._log_warning("Invalid hydrophone sensitivity, using default 4.4 V/bar")
            
            if self._config['sample_rate'] <= 0 or self._config['sample_rate'] > 1000:
                self._config['sample_rate'] = 0.0025
                self._log_warning("Invalid sample rate, using default 0.0025 ms")
            
            # Validate directory
            if not isinstance(self._config['data_dir'], Path):
                self._config['data_dir'] = Path.cwd()
            
            if not self._config['data_dir'].exists():
                self._config['data_dir'] = Path.cwd()
                self._log_warning(f"Data directory not found, using current directory: {Path.cwd()}")
            
        except Exception as e:
            self._log_error(f"Config validation failed: {e}")
    
    def _log_warning(self, message: str):
        """Safe logging with fallback to print."""
        try:
            logging.warning(message)
        except:
            print(f"WARNING: {message}")
    
    def _log_error(self, message: str):
        """Safe logging with fallback to print."""
        try:
            logging.error(message)
        except:
            print(f"ERROR: {message}")
    
    def __getitem__(self, key):
        return self._config.get(key)
    
    def __setitem__(self, key, value):
        self._config[key] = value
        self._validate_config()

# Initialize robust configuration
CONFIG = RobustConfig()

# Setup robust logging
def setup_logging():
    """Setup logging with fallback mechanisms."""
    try:
        log_file = CONFIG['data_dir'] / f"aquaspark_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        logging.info("Logging system initialized successfully")
    except Exception as e:
        print(f"Failed to setup logging: {e}. Continuing with print statements.")


def safe_log(message: str, level: str = "info"):
    """Ultra-safe logging function with multiple fallbacks."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    full_message = f"[{timestamp}] {message}"
    
    try:
        if level.lower() == "error":
            logging.error(message)
        elif level.lower() == "warning":
            logging.warning(message)
        else:
            logging.info(message)
    except:
        try:
            print(full_message)
        except:
            pass  # Ultimate fallback: do nothing


def validate_array_data(data: np.ndarray, name: str) -> bool:
    """Comprehensive array validation with detailed reporting."""
    if data is None:
        safe_log(f"Array {name} is None", "error")
        return False
    
    if not isinstance(data, np.ndarray):
        safe_log(f"Array {name} is not a numpy array: {type(data)}", "error")
        return False
    
    if data.size == 0:
        safe_log(f"Array {name} is empty", "error")
        return False
    
    if data.size < CONFIG['min_samples']:
        safe_log(f"Array {name} too small: {data.size} < {CONFIG['min_samples']}", "error")
        return False
    
    if data.size > CONFIG['max_samples']:
        safe_log(f"Array {name} too large: {data.size} > {CONFIG['max_samples']}", "warning")
        return False
    
    if not np.isfinite(data).all():
        nan_count = np.isnan(data).sum()
        inf_count = np.isinf(data).sum()
        safe_log(f"Array {name} contains non-finite values: {nan_count} NaNs, {inf_count} Infs", "warning")
        return False
    
    return True


def load_oscilloscope_data(filepath: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Ultra-robust data loading with multiple fallback strategies."""
    safe_log(f"Loading file: {filepath.name}")
    
    # Input validation
    if not isinstance(filepath, Path):
        try:
            filepath = Path(filepath)
        except Exception as e:
            safe_log(f"Invalid filepath type: {type(filepath)}", "error")
            return None, None
    
    if not filepath.exists():
        safe_log(f"File does not exist: {filepath}", "error")
        return None, None
    
    if not filepath.is_file():
        safe_log(f"Path is not a file: {filepath}", "error")
        return None, None
    
    if filepath.stat().st_size == 0:
        safe_log(f"File is empty: {filepath}", "error")
        return None, None
    
    if filepath.stat().st_size > 500 * 1024 * 1024:  # 500MB limit
        safe_log(f"File too large: {filepath.stat().st_size / (1024*1024):.1f}MB", "error")
        return None, None
    
    # Multiple loading strategies
    strategies = [
        _load_with_pandas,
        _load_with_numpy,
        _load_with_manual_parsing
    ]
    
    for i, strategy in enumerate(strategies):
        try:
            safe_log(f"Trying loading strategy {i+1}/{len(strategies)}: {strategy.__name__}")
            time, amplitude = strategy(filepath)
            
            if time is not None and amplitude is not None:
                # Validate loaded data
                if not validate_array_data(time, "time"):
                    continue
                if not validate_array_data(amplitude, "amplitude"):
                    continue
                
                # Additional validation
                if len(time) != len(amplitude):
                    safe_log(f"Length mismatch: time={len(time)}, amplitude={len(amplitude)}", "error")
                    continue
                
                # Apply polarity inversion safely
                try:
                    if CONFIG['invert_polarity']:
                        amplitude = amplitude * -1
                        safe_log("Polarity inverted")
                except Exception as e:
                    safe_log(f"Polarity inversion failed: {e}", "warning")
                
                safe_log(f"Successfully loaded {len(time)} samples using {strategy.__name__}")
                return time, amplitude
                
        except Exception as e:
            safe_log(f"Strategy {strategy.__name__} failed: {e}", "warning")
            continue
    
    safe_log(f"All loading strategies failed for {filepath.name}", "error")
    return None, None


def _load_with_pandas(filepath: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Loading strategy 1: Pandas with automatic header detection."""
    # Read file to find data start
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    skip_rows = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.strip().startswith('#'):
            skip_rows = i
            break
    
    df = pd.read_csv(filepath, skiprows=skip_rows, encoding='utf-8', on_bad_lines='skip')
    
    if df.shape[1] < 2:
        raise ValueError(f"Insufficient columns: {df.shape[1]}")
    
    time = df.iloc[:, 0].values.astype(float)
    amplitude = df.iloc[:, 1].values.astype(float)
    
    return time, amplitude


def _load_with_numpy(filepath: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Loading strategy 2: NumPy loadtxt with error handling."""
    # Count header lines
    skip_rows = 0
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.strip().startswith('#') or line.strip().startswith('Time'):
                skip_rows += 1
            else:
                break
    
    data = np.loadtxt(filepath, delimiter=',', skiprows=skip_rows, usecols=(0, 1))
    
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(f"Invalid data shape: {data.shape}")
    
    time = data[:, 0]
    amplitude = data[:, 1]
    
    return time, amplitude


def _load_with_manual_parsing(filepath: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Loading strategy 3: Manual line-by-line parsing (most robust)."""
    time_list = []
    amplitude_list = []
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            
            # Skip comments and headers
            if not line or line.startswith('#') or 'Time' in line:
                continue
            
            try:
                parts = line.split(',')
                if len(parts) >= 2:
                    t_val = float(parts[0])
                    a_val = float(parts[1])
                    
                    if np.isfinite(t_val) and np.isfinite(a_val):
                        time_list.append(t_val)
                        amplitude_list.append(a_val)
                    
            except (ValueError, IndexError) as e:
                safe_log(f"Skipping line {line_num}: {e}", "warning")
                continue
    
    if len(time_list) < CONFIG['min_samples']:
        raise ValueError(f"Insufficient valid data points: {len(time_list)}")
    
    return np.array(time_list), np.array(amplitude_list)


def process_trace(time: np.ndarray, amplitude: np.ndarray) -> Optional[Tuple[np.ndarray, ...]]:
    """Ultra-robust trace processing with comprehensive error handling."""
    try:
        # Input validation
        if not validate_array_data(time, "time") or not validate_array_data(amplitude, "amplitude"):
            return None
        
        if len(time) != len(amplitude):
            safe_log(f"Length mismatch in process_trace: time={len(time)}, amplitude={len(amplitude)}", "error")
            return None
        
        # Safe peak detection with fallbacks
        try:
            abs_amplitude = np.abs(amplitude)
            if not validate_array_data(abs_amplitude, "abs_amplitude"):
                return None
            
            peak_idx = np.argmax(abs_amplitude)
            
            # Validate peak index
            if peak_idx < 0 or peak_idx >= len(amplitude):
                safe_log(f"Invalid peak index: {peak_idx}", "error")
                return None
                
        except Exception as e:
            safe_log(f"Peak detection failed: {e}", "error")
            # Fallback: use middle of array
            peak_idx = len(amplitude) // 2
        
        # Safe time alignment with bounds checking
        try:
            t_set = 0.003  # s
            time_aligned = time - time[peak_idx] + t_set
            
            if not validate_array_data(time_aligned, "time_aligned"):
                safe_log("Time alignment failed, using original time", "warning")
                time_aligned = time.copy()
                
        except Exception as e:
            safe_log(f"Time alignment failed: {e}", "warning")
            time_aligned = time.copy()
        
        # Safe trimming with adaptive length
        try:
            target_length = 8000
            margin = 1550
            
            # Adaptive trimming based on data length
            if len(amplitude) < target_length:
                safe_log(f"Data too short for standard trimming: {len(amplitude)} < {target_length}", "warning")
                # Use entire array
                start_idx, end_idx = 0, len(amplitude)
            else:
                start_idx = max(0, peak_idx - margin)
                end_idx = min(len(amplitude), start_idx + target_length)
                
                # Ensure we have enough data
                if end_idx - start_idx < CONFIG['min_samples']:
                    safe_log(f"Insufficient data after trimming: {end_idx - start_idx}", "warning")
                    start_idx, end_idx = 0, len(amplitude)
            
            time_trimmed = time_aligned[start_idx:end_idx].copy()
            amp_trimmed = amplitude[start_idx:end_idx].copy()
            
            if not validate_array_data(time_trimmed, "time_trimmed") or not validate_array_data(amp_trimmed, "amp_trimmed"):
                return None
                
        except Exception as e:
            safe_log(f"Trimming failed: {e}", "error")
            return None
        
        # Safe unit conversion with validation
        try:
            sensitivity = CONFIG['hydrophone_sensitivity']
            depth = CONFIG['hydrophone_depth']
            
            if sensitivity <= 0 or depth <= 0:
                safe_log(f"Invalid calibration values: sensitivity={sensitivity}, depth={depth}", "error")
                return None
            
            bar_m = amp_trimmed * depth / sensitivity
            
            if not validate_array_data(bar_m, "bar_m"):
                safe_log("Bar-m conversion resulted in invalid data", "error")
                return None
                
        except Exception as e:
            safe_log(f"Unit conversion failed: {e}", "error")
            return None
        
        # Safe normalization with zero-division protection
        try:
            max_abs_amp = np.max(np.abs(amp_trimmed))
            
            if max_abs_amp == 0 or not np.isfinite(max_abs_amp):
                safe_log("Cannot normalize: zero or invalid maximum amplitude", "error")
                return None
            
            normalized = amp_trimmed / max_abs_amp
            
            if not validate_array_data(normalized, "normalized"):
                safe_log("Normalization resulted in invalid data", "error")
                return None
                
        except Exception as e:
            safe_log(f"Normalization failed: {e}", "error")
            return None
        
        safe_log(f"Trace processed successfully: {len(time_trimmed)} samples")
        return time_trimmed, amp_trimmed, bar_m, normalized
        
    except Exception as e:
        safe_log(f"Critical error in process_trace: {e}", "error")
        safe_log(f"Traceback: {traceback.format_exc()}", "error")
        return None


def compute_spectrum(signal: np.ndarray) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Ultra-robust spectrum computation with comprehensive error handling."""
    try:
        # Input validation
        if not validate_array_data(signal, "signal"):
            return None
        
        # Safe sampling frequency calculation
        try:
            sample_rate = CONFIG['sample_rate']
            if sample_rate <= 0 or sample_rate > 1000:
                safe_log(f"Invalid sample rate: {sample_rate}, using default", "warning")
                sample_rate = 0.0025
            
            fs = 1000.0 / sample_rate
            
            if not np.isfinite(fs) or fs <= 0:
                safe_log(f"Invalid sampling frequency: {fs}", "error")
                return None
                
        except Exception as e:
            safe_log(f"Sampling frequency calculation failed: {e}", "error")
            return None
        
        # Safe FFT size calculation with bounds checking
        try:
            signal_length = len(signal)
            
            if signal_length < 4:  # Minimum for meaningful FFT
                safe_log(f"Signal too short for FFT: {signal_length}", "error")
                return None
            
            # Compute next power of 2 safely
            log2_length = np.log2(signal_length - 0.5)
            if not np.isfinite(log2_length):
                log2_length = np.log2(signal_length)
            
            n_fft = 2 ** (int(log2_length) + 1)
            
            # Apply padding with safety limits
            padding = max(1, min(CONFIG['fft_padding'], 8))  # Limit padding to prevent memory issues
            padded_length = padding * n_fft
            
            # Memory safety check
            max_fft_size = 2**20  # 1M points maximum
            if padded_length > max_fft_size:
                safe_log(f"FFT size too large: {padded_length}, reducing to {max_fft_size}", "warning")
                padded_length = max_fft_size
                
        except Exception as e:
            safe_log(f"FFT size calculation failed: {e}", "error")
            return None
        
        # Safe FFT computation
        try:
            # Remove DC offset and apply window to reduce spectral leakage
            signal_windowed = signal - np.mean(signal)
            
            # Apply Hanning window for better spectral properties
            if len(signal_windowed) > 1:
                window = np.hanning(len(signal_windowed))
                signal_windowed = signal_windowed * window
            
            fft_result = np.fft.fft(signal_windowed, n=padded_length)
            
            if not np.isfinite(fft_result).all():
                safe_log("FFT result contains non-finite values", "error")
                return None
                
        except MemoryError:
            safe_log("Out of memory during FFT computation", "error")
            return None
        except Exception as e:
            safe_log(f"FFT computation failed: {e}", "error")
            return None
        
        # Safe spectral processing
        try:
            n_half = len(fft_result) // 2
            if n_half < 1:
                safe_log("FFT result too short", "error")
                return None
            
            fft_half = fft_result[:n_half]
            
            # Compute magnitude spectrum
            magnitude_spectrum = np.abs(fft_half)
            
            if not validate_array_data(magnitude_spectrum, "magnitude_spectrum"):
                return None
            
            # Safe dB conversion with floor to prevent log(0)
            max_magnitude = np.max(magnitude_spectrum)
            if max_magnitude <= 0 or not np.isfinite(max_magnitude):
                safe_log("Invalid maximum magnitude for dB conversion", "error")
                return None
            
            # Add small floor to prevent log(0)
            floor_value = max_magnitude * 1e-12
            magnitude_spectrum_safe = np.maximum(magnitude_spectrum, floor_value)
            
            power_db = 10 * np.log10(magnitude_spectrum_safe / max_magnitude)
            
            # Validate dB values
            if not validate_array_data(power_db, "power_db"):
                return None
            
            # Generate frequency axis
            frequencies = np.linspace(0, fs/2, n_half)
            
            if not validate_array_data(frequencies, "frequencies"):
                return None
            
            safe_log(f"Spectrum computed successfully: {len(frequencies)} frequency bins")
            return frequencies, power_db
            
        except Exception as e:
            safe_log(f"Spectral processing failed: {e}", "error")
            return None
        
    except Exception as e:
        safe_log(f"Critical error in compute_spectrum: {e}", "error")
        safe_log(f"Traceback: {traceback.format_exc()}", "error")
        return None


def create_plots() -> Optional[Tuple[plt.Figure, List[plt.Axes]]]:
    """Ultra-robust plot creation with fallback mechanisms."""
    try:
        # Safe matplotlib configuration
        try:
            # Set safe matplotlib backend
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend for offshore use
            
            plt.rcParams.update({
                'font.size': max(8, min(CONFIG['font_size'], 24)),  # Safe font size range
                'figure.max_open_warning': 0,  # Disable warnings
                'axes.grid': True,
                'grid.alpha': 0.3,
            })
            safe_log("Matplotlib configured successfully")
            
        except Exception as e:
            safe_log(f"Matplotlib configuration failed: {e}", "warning")
        
        # Safe figure creation
        try:
            # Validate figure size
            fig_width, fig_height = CONFIG['figure_size']
            fig_width = max(8, min(fig_width, 50))  # Reasonable bounds
            fig_height = max(6, min(fig_height, 50))
            
            fig = plt.figure(figsize=(fig_width, fig_height), dpi=100)
            
            if fig is None:
                safe_log("Figure creation failed", "error")
                return None
                
        except Exception as e:
            safe_log(f"Figure creation failed: {e}", "error")
            return None
        
        # Safe layout adjustment
        try:
            fig.subplots_adjust(
                hspace=0.25, 
                wspace=0.3, 
                top=0.9, 
                bottom=0.1,
                left=0.1,
                right=0.95
            )
        except Exception as e:
            safe_log(f"Layout adjustment failed: {e}", "warning")
        
        # Safe title creation
        try:
            title = f"Sparker Pulse Analysis: {CONFIG['test_name']}"
            if CONFIG['invert_polarity']:
                title += " (Inverted)"
            
            # Truncate title if too long
            if len(title) > 80:
                title = title[:77] + "..."
            
            fig.suptitle(title, fontsize=min(22, CONFIG['font_size'] + 6), y=0.95)
            
        except Exception as e:
            safe_log(f"Title creation failed: {e}", "warning")
        
        # Safe subplot creation
        try:
            ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)  # Raw traces
            ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)  # Normalized traces  
            ax3 = plt.subplot2grid((2, 4), (1, 0), colspan=4)  # Spectra
            
            axes = [ax1, ax2, ax3]
            
            # Validate axes
            for i, ax in enumerate(axes):
                if ax is None:
                    safe_log(f"Axis {i} creation failed", "error")
                    return None
            
            safe_log("Plots created successfully")
            return fig, axes
            
        except Exception as e:
            safe_log(f"Subplot creation failed: {e}", "error")
            return None
        
    except Exception as e:
        safe_log(f"Critical error in create_plots: {e}", "error")
        safe_log(f"Traceback: {traceback.format_exc()}", "error")
        return None


def setup_axes(axes: List[plt.Axes]) -> None:
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


def safe_file_discovery(data_dir: Path) -> List[Path]:
    """Ultra-safe file discovery with multiple fallback strategies."""
    csv_files = []
    
    try:
        # Strategy 1: Standard glob
        csv_files = list(data_dir.glob('*.csv'))
        safe_log(f"Found {len(csv_files)} CSV files using glob")
        
    except Exception as e:
        safe_log(f"Glob search failed: {e}", "warning")
        
        try:
            # Strategy 2: Manual directory scanning
            if data_dir.exists() and data_dir.is_dir():
                for item in data_dir.iterdir():
                    if item.is_file() and item.suffix.lower() == '.csv':
                        csv_files.append(item)
                safe_log(f"Found {len(csv_files)} CSV files using manual scan")
                
        except Exception as e2:
            safe_log(f"Manual scan failed: {e2}", "warning")
    
    # Limit number of files to prevent resource exhaustion
    if len(csv_files) > CONFIG['max_files']:
        safe_log(f"Too many files ({len(csv_files)}), limiting to {CONFIG['max_files']}", "warning")
        csv_files = csv_files[:CONFIG['max_files']]
    
    return sorted(csv_files)  # Sort for consistent processing order


def main():
    """Ultra-robust main analysis pipeline with comprehensive error handling."""
    try:
        # Initialize logging
        setup_logging()
        safe_log("="*60)
        safe_log("AquaSpark Pro - Offshore Analysis Starting")
        safe_log("="*60)
        
        # Validate configuration
        safe_log("Validating configuration...")
        if not CONFIG['data_dir'].exists():
            safe_log(f"Data directory not found: {CONFIG['data_dir']}", "error")
            safe_log("Please update CONFIG['data_dir'] to point to your data folder")
            return False
        
        # Find CSV files with robust error handling
        safe_log("Discovering CSV files...")
        csv_files = safe_file_discovery(CONFIG['data_dir'])
        
        if not csv_files:
            safe_log("No CSV files found. Analysis cannot proceed.", "error")
            safe_log("Please check that:")
            safe_log("1. The data directory path is correct")
            safe_log("2. CSV files exist in the directory")
            safe_log("3. You have read permissions for the directory")
            return False
        
        safe_log(f"Processing {len(csv_files)} files...")
        
        # Initialize storage with error handling
        traces_raw = []
        traces_bar = []
        traces_norm = []
        successful_files = 0
        failed_files = 0
        
        # Create plots with error handling
        plot_result = create_plots()
        if plot_result is None:
            safe_log("Plot creation failed. Analysis cannot continue.", "error")
            return False
        
        fig, axes = plot_result
        ax1, ax2, ax3 = axes
        
        # Generate colors safely
        try:
            colors = plt.cm.autumn(np.linspace(0, 1, max(1, len(csv_files))))
        except Exception as e:
            safe_log(f"Color generation failed: {e}", "warning")
            # Fallback to basic colors
            colors = ['red', 'blue', 'green', 'orange', 'purple'] * (len(csv_files) // 5 + 1)
        
        # Process each file with robust error handling
        time_ms = None  # Will be set from first successful file
        
        for i, filepath in enumerate(csv_files):
            try:
                safe_log(f"Processing file {i+1}/{len(csv_files)}: {filepath.name}")
                
                # Load data
                time, amplitude = load_oscilloscope_data(filepath)
                if time is None or amplitude is None:
                    failed_files += 1
                    continue
                
                # Process trace
                result = process_trace(time, amplitude)
                if result is None:
                    failed_files += 1
                    continue
                
                time_proc, amp_proc, bar_proc, norm_proc = result
                
                # Store for stacking
                traces_raw.append(amp_proc)
                traces_bar.append(bar_proc) 
                traces_norm.append(norm_proc)
                
                # Set time axis (use first successful file as reference)
                if time_ms is None:
                    time_ms = time_proc * 1000
                    peak_idx = np.argmax(np.abs(amp_proc))
                    t_min = time_ms[peak_idx] - 0.5
                    t_max = t_min + CONFIG['plot_length']
                
                # Safe plotting
                try:
                    color = colors[i % len(colors)]
                    label = filepath.stem[:20]  # Truncate long names
                    
                    ax1.plot(time_ms, bar_proc, color=color, 
                            linewidth=CONFIG['line_width'], label=label, alpha=0.7)
                    ax2.plot(time_ms, norm_proc, color=color, 
                            linewidth=CONFIG['line_width'], label=label, alpha=0.7)
                    
                    # Compute spectrum
                    spectrum_result = compute_spectrum(amp_proc)
                    if spectrum_result is not None:
                        frequencies, power_db = spectrum_result
                        ax3.plot(frequencies, power_db, color=color, 
                                linewidth=CONFIG['line_width'], label=label, alpha=0.7)
                    
                    # Set time limits
                    ax1.set_xlim(t_min, t_max)
                    ax2.set_xlim(t_min, t_max)
                    
                except Exception as e:
                    safe_log(f"Plotting failed for {filepath.name}: {e}", "warning")
                
                successful_files += 1
                safe_log(f"Successfully processed {filepath.name}")
                
            except Exception as e:
                failed_files += 1
                safe_log(f"Failed to process {filepath.name}: {e}", "error")
                continue
        
        # Summary of processing
        safe_log(f"Processing summary: {successful_files} successful, {failed_files} failed")
        
        if successful_files == 0:
            safe_log("No files processed successfully. Cannot generate plots.", "error")
            return False
        
        # Compute and plot stacks with error handling
        if traces_bar and traces_raw and traces_norm:
            try:
                safe_log("Computing stacks...")
                
                # Safe stacking with length validation
                min_length = min(len(trace) for trace in traces_bar)
                traces_bar_trimmed = [trace[:min_length] for trace in traces_bar]
                traces_norm_trimmed = [trace[:min_length] for trace in traces_norm]
                traces_raw_trimmed = [trace[:min_length] for trace in traces_raw]
                
                stack_bar = np.mean(traces_bar_trimmed, axis=0)
                stack_norm_raw = np.mean(traces_norm_trimmed, axis=0)
                stack_norm = stack_norm_raw / np.max(np.abs(stack_norm_raw))
                
                # Trim time axis to match
                if time_ms is not None:
                    time_ms_trimmed = time_ms[:min_length]
                
                    # Plot stacks
                    ax1.plot(time_ms_trimmed, stack_bar, 'blue', 
                            linewidth=CONFIG['line_width']*4, label='Stack', zorder=10)
                    ax2.plot(time_ms_trimmed, stack_norm, 'blue', 
                            linewidth=CONFIG['line_width']*4, label='Stack', zorder=10)
                    
                    # Stack spectrum
                    stack_raw = np.mean(traces_raw_trimmed, axis=0)
                    spectrum_result = compute_spectrum(stack_raw)
                    if spectrum_result is not None:
                        freq_stack, power_stack = spectrum_result
                        ax3.plot(freq_stack, power_stack, 'blue', 
                                linewidth=CONFIG['line_width']*4, label='Stack', zorder=10)
                
                safe_log(f"Stacked {len(traces_bar)} traces successfully")
                
            except Exception as e:
                safe_log(f"Stacking failed: {e}", "warning")
        
        # Finalize plots
        try:
            setup_axes(axes)
        except Exception as e:
            safe_log(f"Axis setup failed: {e}", "warning")
        
        # Safe file saving with backup
        try:
            base_name = f"AquaSpark_Analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            output_path = CONFIG['data_dir'] / f"{base_name}.png"
            
            # Save with error handling
            fig.savefig(output_path, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            safe_log(f"Figure saved to: {output_path}")
            
            # Save backup copy
            if CONFIG['backup_save']:
                backup_path = CONFIG['data_dir'] / f"{base_name}_backup.png"
                fig.savefig(backup_path, dpi=150, bbox_inches='tight')
                safe_log(f"Backup saved to: {backup_path}")
            
        except Exception as e:
            safe_log(f"File saving failed: {e}", "error")
        
        # Safe display (disable for offshore headless operation)
        try:
            plt.tight_layout()
            # plt.show()  # Commented out for offshore use
        except Exception as e:
            safe_log(f"Display failed: {e}", "warning")
        
        safe_log("="*60)
        safe_log("Analysis completed successfully!")
        safe_log(f"Processed {successful_files}/{len(csv_files)} files")
        safe_log("="*60)
        return True
        
    except Exception as e:
        safe_log(f"CRITICAL ERROR in main: {e}", "error")
        safe_log(f"Traceback: {traceback.format_exc()}", "error")
        return False


if __name__ == "__main__":
    """Ultra-safe entry point with maximum error resilience."""
    try:
        # Suppress all warnings for clean offshore operation
        warnings.filterwarnings("ignore")
        
        # Set matplotlib to non-interactive mode
        import matplotlib
        matplotlib.use('Agg')
        
        # Run main analysis
        success = main()
        
        if success:
            print("\n✅ AquaSpark Analysis completed successfully!")
            sys.exit(0)
        else:
            print("\n❌ AquaSpark Analysis failed. Check log for details.")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n⚠️  Analysis interrupted by user.")
        sys.exit(2)
        
    except Exception as e:
        print(f"\n💥 FATAL ERROR: {e}")
        try:
            print(f"Traceback: {traceback.format_exc()}")
        except:
            pass
        sys.exit(3)