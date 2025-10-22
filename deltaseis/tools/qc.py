"""
Seismic Data Quality Control (QC) Module

This module provides functions for analyzing the quality of seismic data using
various metrics including signal-to-noise ratios and wavelet entropy measures.

Main Functions:
--------------
- compute_windowed_snr: Time-domain SNR calculation
- compute_frequency_snr: Frequency-domain SNR calculation  
- compute_wavelet_entropy: Wavelet-based complexity measure
- analyze_seismic_data: Complete QC analysis with PDF report generation

Example:
--------
>>> from deltaseis.tools.qc import analyze_seismic_data
>>> import numpy as np
>>> 
>>> # Load seismic data (traces x samples)
>>> data = np.random.randn(1000, 2000)  # Example data
>>> fs = 4000  # Sampling frequency in Hz
>>> 
>>> # Define analysis windows
>>> signal_window = (500, 1500)  # Signal of interest
>>> noise_window = (0, 400)      # Background noise
>>> 
>>> # Generate QC report
>>> analyze_seismic_data(data, fs, signal_window, noise_window, 
>>>                     output_pdf='qc_report.pdf', max_traces=100)

Authors: Roeland Nieboer, Deltares
"""

import numpy as np
import matplotlib.pyplot as plt
import pywt
from scipy.signal import welch
from matplotlib.backends.backend_pdf import PdfPages

def compute_windowed_snr(trace, signal_window, noise_window):
    """
    Compute Signal-to-Noise Ratio (SNR) using time-domain windows.
    
    Parameters:
    -----------
    trace : array_like
        Input seismic trace data
    signal_window : tuple
        (start, end) sample indices for signal window
    noise_window : tuple
        (start, end) sample indices for noise window
        
    Returns:
    --------
    float
        SNR in decibels (dB)
        
    Notes:
    ------
    SNR is calculated as 10 * log10(signal_power / noise_power)
    where power is the mean squared amplitude in each window.
    """
    signal = trace[signal_window[0]:signal_window[1]]
    noise = trace[noise_window[0]:noise_window[1]]
    
    signal_power = np.mean(signal**2)
    noise_power = np.mean(noise**2)
    
    # Handle edge cases
    if noise_power == 0 or signal_power == 0:
        return 0.0  # Return 0 dB for edge cases
    
    snr = 10 * np.log10(signal_power / noise_power)
    
    # Clip extreme values to reasonable range
    snr = np.clip(snr, -100, 100)
    
    return snr

def compute_frequency_snr(trace, fs, signal_band=(100, 1000), noise_band=(0, 50)):
    """
    Compute Signal-to-Noise Ratio (SNR) in the frequency domain.
    
    Parameters:
    -----------
    trace : array_like
        Input seismic trace data
    fs : float
        Sampling frequency in Hz
    signal_band : tuple, optional
        (low_freq, high_freq) in Hz for signal band (default: 100-1000 Hz)
    noise_band : tuple, optional
        (low_freq, high_freq) in Hz for noise band (default: 0-50 Hz)
        
    Returns:
    --------
    float
        SNR in decibels (dB)
        
    Notes:
    ------
    Uses Welch's method to estimate power spectral density, then computes
    SNR as the ratio of power in signal band to power in noise band.
    """
    f, Pxx = welch(trace, fs=fs, nperseg=256)
    signal_power = np.sum(Pxx[(f >= signal_band[0]) & (f <= signal_band[1])])
    noise_power = np.sum(Pxx[(f >= noise_band[0]) & (f <= noise_band[1])])
    
    # Handle edge cases
    if noise_power == 0 or signal_power == 0:
        return 0.0  # Return 0 dB for edge cases
    
    snr = 10 * np.log10(signal_power / noise_power)
    
    # Clip extreme values to reasonable range
    snr = np.clip(snr, -100, 100)
    
    return snr

def compute_wavelet_entropy(trace, wavelet='db4', level=4):
    """
    Compute wavelet entropy as a measure of signal complexity.
    
    Parameters:
    -----------
    trace : array_like
        Input seismic trace data
    wavelet : str, optional
        Wavelet type for decomposition (default: 'db4' - Daubechies 4)
    level : int, optional
        Number of decomposition levels (default: 4)
        
    Returns:
    --------
    float
        Wavelet entropy value
        
    Notes:
    ------
    Wavelet entropy measures the degree of order/disorder in the signal.
    Lower entropy indicates more structured/coherent signals, while
    higher entropy indicates more random/noisy signals.
    
    The entropy is calculated as:
    H = -sum(p_i * log2(p_i))
    where p_i is the normalized energy in each wavelet coefficient band.
    """
    coeffs = pywt.wavedec(trace, wavelet, level=level)
    energy = np.array([np.sum(c**2) for c in coeffs])
    total_energy = np.sum(energy)
    prob = energy / total_energy
    entropy = -np.sum(prob * np.log2(prob + 1e-12))  # avoid log(0)
    return entropy

def analyze_seismic_data(matrix, fs, signal_window, noise_window, output_pdf='seismic_quality_report.pdf', max_traces=100):
    """
    Fast seismic QC analysis - samples traces for efficiency
    
    Parameters:
    -----------
    matrix : array_like
        Seismic data matrix (traces x samples)
    fs : float
        Sampling frequency
    signal_window : tuple
        (start, end) samples for signal window
    noise_window : tuple
        (start, end) samples for noise window  
    output_pdf : str
        Output PDF filename
    max_traces : int
        Maximum number of traces to analyze (default 100 for speed)
    """
    
    n_traces, n_samples = matrix.shape
    print(f"Analyzing seismic data: {n_traces} traces x {n_samples} samples")
    
    # Sample traces if file is too large
    if n_traces > max_traces:
        print(f"Sampling {max_traces} traces from {n_traces} total traces...")
        trace_indices = np.linspace(0, n_traces-1, max_traces, dtype=int)
        sampled_matrix = matrix[trace_indices]
    else:
        trace_indices = np.arange(n_traces)
        sampled_matrix = matrix
        
    print("Computing QC metrics...")
    
    # Vectorized computation for speed
    snr_time_list = []
    snr_freq_list = []
    entropy_list = []
    
    for i, trace in enumerate(sampled_matrix):
        if i % 20 == 0:  # Progress indicator
            print(f"Processing trace {i+1}/{len(sampled_matrix)}")
            
        snr_time = compute_windowed_snr(trace, signal_window, noise_window)
        snr_freq = compute_frequency_snr(trace, fs)
        entropy = compute_wavelet_entropy(trace)
        
        snr_time_list.append(snr_time)
        snr_freq_list.append(snr_freq)
        entropy_list.append(entropy)
    
    # Convert to arrays for statistics and filter out invalid values
    snr_time_array = np.array(snr_time_list)
    snr_freq_array = np.array(snr_freq_list)
    entropy_array = np.array(entropy_list)
    
    # Filter out infinite and NaN values
    valid_mask = (np.isfinite(snr_time_array) & 
                  np.isfinite(snr_freq_array) & 
                  np.isfinite(entropy_array))
    
    snr_time_array = snr_time_array[valid_mask]
    snr_freq_array = snr_freq_array[valid_mask]
    entropy_array = entropy_array[valid_mask]
    
    if len(snr_time_array) == 0:
        print("Warning: No valid traces found after filtering. Skipping QC report generation.")
        return
    
    print(f"Generating QC report for {len(snr_time_array)} valid traces...")
    
    with PdfPages(output_pdf) as pdf:
        # Page 1: Summary Statistics
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(f'Seismic QC Summary - {n_traces} traces ({len(sampled_matrix)} analyzed)', fontsize=14)
        
        # SNR Time Domain Histogram
        axs[0,0].hist(snr_time_array, bins=20, alpha=0.7, edgecolor='black')
        axs[0,0].set_title(f'Time SNR Distribution\nMean: {np.mean(snr_time_array):.1f} dB')
        axs[0,0].set_xlabel('SNR (dB)')
        axs[0,0].set_ylabel('Count')
        
        # SNR Frequency Domain Histogram  
        axs[0,1].hist(snr_freq_array, bins=20, alpha=0.7, edgecolor='black', color='orange')
        axs[0,1].set_title(f'Freq SNR Distribution\nMean: {np.mean(snr_freq_array):.1f} dB')
        axs[0,1].set_xlabel('SNR (dB)')
        axs[0,1].set_ylabel('Count')
        
        # Entropy Histogram
        axs[1,0].hist(entropy_array, bins=20, alpha=0.7, edgecolor='black', color='green')
        axs[1,0].set_title(f'Wavelet Entropy Distribution\nMean: {np.mean(entropy_array):.2f}')
        axs[1,0].set_xlabel('Entropy')
        axs[1,0].set_ylabel('Count')
        
        # Trace Quality vs Position (only valid traces)
        valid_trace_indices = trace_indices[valid_mask] if n_traces > len(sampled_matrix) else np.arange(len(snr_time_array))
        axs[1,1].plot(valid_trace_indices, snr_time_array, 'o-', alpha=0.7, markersize=3)
        axs[1,1].set_title('SNR vs Trace Position')
        axs[1,1].set_xlabel('Trace Number')
        axs[1,1].set_ylabel('Time SNR (dB)')
        axs[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        
        # Page 2: Sample Traces (show first, middle, last few traces)
        sample_indices = [0, len(sampled_matrix)//4, len(sampled_matrix)//2, 
                         3*len(sampled_matrix)//4, len(sampled_matrix)-1]
        
        fig, axs = plt.subplots(len(sample_indices), 1, figsize=(12, 10))
        fig.suptitle('Sample Traces', fontsize=14)
        
        time_axis = np.arange(n_samples) / fs * 1000  # Convert to milliseconds
        
        for i, idx in enumerate(sample_indices):
            if idx < len(sampled_matrix):
                trace = sampled_matrix[idx]
                actual_trace_num = trace_indices[idx]
                axs[i].plot(time_axis, trace, 'b-', linewidth=0.5)
                axs[i].set_title(f'Trace {actual_trace_num} - SNR: {snr_time_list[idx]:.1f} dB')
                axs[i].set_ylabel('Amplitude')
                if i == len(sample_indices)-1:
                    axs[i].set_xlabel('Time (ms)')
                axs[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
    
    print(f"QC report saved to: {output_pdf}")
    print(f"Summary: SNR = {np.mean(snr_time_array):.1f}±{np.std(snr_time_array):.1f} dB, "
          f"Entropy = {np.mean(entropy_array):.2f}±{np.std(entropy_array):.2f}")

