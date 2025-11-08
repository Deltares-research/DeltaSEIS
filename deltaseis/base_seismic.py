# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 09:14:16 2022
Seismic Class contaiing tools for seismic processing
Currently contains:
- agc gain,
- bandpass_filter 
- f-k spectrum, 
- plot f-k spectrum
- plot histogram
- plot_psd

UNDER CONSTRUCTION:
- plot_time_psd_fk OR can we cast results from other plot funcitons in a pyplot 3 axs??
---- return plots as im? or fig?

to be added:
- cut record length
- f-k filter, use bowpy?
- interactive line drawing for slopes
- xarray?
- fbpick for shot selection
- add wiggle module (using segpy)
- add decon module (using thrid party)?
- add migration module (using thrid party)?
- add multiple suppression module (using thrid party)?
- add AVO module (using thrid party)?
- add inversion module (using thrid party)?
- add QC module (using thrid party)?
- add kwargs??
- t2 gain

- nmo (stretch mute?)
- stack

Readers: TDMS, sg2, dat, segy, seg, segd have to be made in separately
Adding tests becomes more important

@author: Roeland Nieboer @ Deltares
"""

import matplotlib.pyplot as plt
from matplotlib import colors 
import numpy as np
from scipy import signal
from scipy.fft import fft2, fftshift
from scipy.ndimage import uniform_filter1d
import sys
          
class Seismic:
    '''Tools for basic processing and visualization of seismic data.'''
      
    def __init__(self, data, fs, dx):
             
        self.data = data
        self.fs = fs
        self.dx = dx
        self.byte_format = str(data.dtype)
        self.time_vector = np.linspace(0, self.data.shape[0]/self.fs, self.data.shape[0])

    def convert_to_trace_data(self, data_sample_format=None, inplace=True):
        """
        Convert 2D array to a list of trace arrays to be written in SEGY format.
        This method converts the internal 2D data array into a list of trace arrays,
        which can then be written in SEGY format. The data sampling format can be 
        specified to ensure that no narrowing of the data occurs when writing to SEGY.
        Parameters:
        -----------
        data_sample_format : str, optional
            The desired data sample format. Supported formats are 'int16', 'int32', 
            'float32', and 'float64'. If None, the data will not be converted.
        inplace : bool, optional
            If True, the converted data will replace the original data in the object.
            If False, the converted data will be returned. Default is True.
        Returns:
        --------
        list of numpy.ndarray or None
            If inplace is False, returns a list of trace arrays. If inplace is True, 
            the method returns None and updates the object's data attribute.
        Raises:
        -------
        ValueError
            If an unsupported data_sample_format is provided.
        """

        if data_sample_format is not None:
            formats = {
                'int16': np.int16,
                'int32': np.int32,
                'float32': np.float32,
                'float64': np.float64
            }
            if data_sample_format in formats:
                if 'int' in data_sample_format:
                    data = (self.data * np.iinfo(formats[data_sample_format]).max / np.max(np.abs(self.data))).astype(formats[data_sample_format])
                else:
                    data = self.data.astype(formats[data_sample_format])
            else:
                raise ValueError("Unsupported data_sample_format. Use 'int16', 'int32', 'float32', or 'float64'.")
        else:
            data = self.data

        data = [np.array(data.T[trace]) for trace in range(data.shape[1])]

        if inplace:
            self.data = data
        else:
            return data
    
    def agc_gain(self, data=None, agc_gate=0.25, method='rms', inplace=True):
        
        '''
        Applying AGC gain to seismic data. Fast way to find trends in the F-K domain
        However, it is not the same as true amplitude correction, so the input gives 
        no valid results when using for AVO analysis or inversion.
              
        Parameters
        ----------
        data: 2d array, optional
            timeseries of seismic data, default=None
            if not defined uses the data variable of the instance
        agc_gate : float, optional
            length of window in seconds to use for the AGC gain correction.
            The default is 0.25 s
        method : string, optional
            Normalization method for AGC correction. The default is 'rms'.
    
        Returns
        -------
        data_agc : 2D array of floats
            seismic data that has AGC gain applied. 
    
        '''
        
        if data is None:
            self.data64 = self.data.astype(np.float64)    
        else:
            self.data64 = data.astype(np.float64)
            
        samples, traces = np.shape(self.data64)
        
        samples_agc = (agc_gate*self.fs) + 1
        samples_agc = np.floor(samples_agc/2)
        
        triangle = signal.windows.triang(2*samples_agc+1)
        
        data_agc = np.zeros((np.shape(self.data64)))
        
        trace_squared = np.power(self.data64, 2)
        rms = np.sqrt(signal.convolve2d(trace_squared, triangle[:, np.newaxis], mode='same'))
        epsi = 1.e-10 * np.max(rms, axis=0)
        op = rms / (np.power(rms, 2) + epsi)
        data_agc = self.data64 * op
        
        if method=='amplitude':            # Normalize by amplitude
        
            max_amplitude = np.max(np.abs(data_agc), axis=0)
            data_agc /= max_amplitude
        
        elif method=='rms':                # Normalize by rms
            
            rms_amplitude = np.sqrt(np.sum(data_agc**2, axis=0) / samples)
            data_agc /= rms_amplitude
        
        else:
            print("please specify method as 'rms'(default) or 'amplitude'")
            
        if inplace==True:
            self.data = data_agc
        else:
            return data_agc 
    
    def time_power_gain(self, power=2, inplace=True, data=None):
        """
        Apply time power (T^power) gain to the seismic data.
        This method applies a time power gain to the seismic data, which can enhance
        the visibility of deeper reflections in the data.
        inplace : bool, optional
            If True, modify the instance's data in place. If False, return a new array
            with the time power gain applied. Default is True.
        data : 2D array, optional
            Timeseries of seismic data. If not provided, the instance's data attribute
            will be used. Default is None.
            Seismic data with time power gain applied. Only returned if inplace is False.
        Notes
        -----
        The time power gain is applied by multiplying each data point by the power of
        its corresponding time value.
        """
    
        if data is None:
            data = self.data

        data_t2_gain = (self.time_vector[:, np.newaxis] ** power) * data

        if inplace==True:
            self.data = data_t2_gain
        else:       
            return data_t2_gain
        
    def spherical_divergence_gain(self, rms_velocity, inplace=True, data=None):
        '''
        This correction compensates for loss of amplitudes due to
        spherical wavefront spreading.

        Parameters
        ----------
        rms_velocity : float or 1D array
            RMS velocity. If a scalar is provided, it is assumed constant for all times.
            If a vector is provided, it should have the same length as the number of samples.
        inplace : bool, optional
            If True, modify the instance's data in place. If False, return a new array
            with the spherical divergence gain applied. Default is True.
        data : 2D array, optional
            Timeseries of seismic data. If not provided, the instance's data attribute
            will be used. Default is None.

        Returns
        -------
        data_spherical_divergence_gain : 2D array
            Seismic data with spherical divergence gain applied. Only returned if inplace is False.
        '''
        if data is None:
            data = self.data

        if np.isscalar(rms_velocity):
            rms_velocity = np.full(len(self.time_vector), rms_velocity)

        data_spherical_divergence_gain = (self.time_vector[:, np.newaxis] * rms_velocity[:, np.newaxis]) * data

        if inplace:
            self.data = data_spherical_divergence_gain
        else:
            return data_spherical_divergence_gain

    def cylindrical_divergence_gain(self, rms_velocity, inplace=True, data=None):
        '''
        This correction compensates for loss of amplitudes due to
        cylindrical wavefront spreading.

        Parameters
        ----------
        rms_velocity : float or 1D array
            RMS velocity. If a scalar is provided, it is assumed constant for all times.
            If a vector is provided, it should have the same length as the number of samples.
        inplace : bool, optional
            If True, modify the instance's data in place. If False, return a new array
            with the spherical divergence gain applied. Default is True.
        data : 2D array, optional
            Timeseries of seismic data. If not provided, the instance's data attribute
            will be used. Default is None.

        Returns
        -------
        data_cylindrical_divergence_gain : 2D array
            Seismic data with spherical divergence gain applied. Only returned if inplace is False.
        '''
        if data is None:
            data = self.data

        if np.isscalar(rms_velocity):
            rms_velocity = np.full(len(self.time_vector), rms_velocity)

        data_cylindrical_divergence_gain = (self.time_vector[:, np.newaxis] * rms_velocity[:, np.newaxis]) * data

        if inplace:
            self.data = data_cylindrical_divergence_gain
        else:
            return data_cylindrical_divergence_gain
        
    def anelastic_attenuation_gain(self, anelestic_attenuation_factor, rms_velocity, inplace=True, data=None):
        '''
        This correction compensates for loss of amplitudes due to anelastic attenuation. This is due
        to fluid movement in th pores of rock (sloshing) and grain boundary friction (jostling).
        Sloshing has a stronger effect on attenuation than jostling.

        Parameters
        ----------
        rms_velocity : float or 1D array
            RMS velocity. If a scalar is provided, it is assumed constant for all times.
            If a vector is provided, it should have the same length as the number of samples.
        inplace : bool, optional
            If True, modify the instance's data in place. If False, return a new array
            with the spherical divergence gain applied. Default is True.
        data : 2D array, optional
            Timeseries of seismic data. If not provided, the instance's data attribute
            will be used. Default is None.

        Returns
        -------
        data_anelastic_attenuation_gain : 2D array
            Seismic data with spherical divergence gain applied. Only returned if inplace is False.
        '''
        if data is None:
            data = self.data

        print(f"An anelastic attenuation factor of {anelestic_attenuation_factor} corresponds to a Q value of {1 / anelestic_attenuation_factor}")
        self.anelestic_attenuation_factor = anelestic_attenuation_factor
        self.q = 1 / anelestic_attenuation_factor

        if np.isscalar(rms_velocity):
            rms_velocity = np.full(len(self.time_vector), rms_velocity)

        data_anelastic_attenuation_gain = np.exp(-anelestic_attenuation_factor * self.time_vector[:, np.newaxis] * rms_velocity[:, np.newaxis]) * data

        if inplace:
            self.data = data_anelastic_attenuation_gain
        else:
            return data_anelastic_attenuation_gain
    
    def trace_averaging(self, number_of_traces=1, inplace=True):
        '''
        Calculates a moving average over 2 x n +1 trace, so n traces on both sides
        of each trace. At the start trace(s) the end trace(s) are added, idem for 
        end trace(s) which include start trace(s). 
        
        Parameters
        ----------
        number_of_traces : int, optional
            Number of 'neighboring' traces on both sides to use in averaging

        Returns
        -------
        data : 2D Array of floats
            2D trace-averaged data
        '''
        data = uniform_filter1d(self.data, size=2 * number_of_traces + 1, axis=1, mode='wrap')
        
        if inplace==True:
            self.data = data
        else:
            return data
        
    def bandpass_filter(self, lowcut, highcut, forder=3, inplace=True):
        
        """
            Applies a bandpass filter to the inputted signal, leaving only the signal's frequency content between the highcut and lowcut frequencies and filtering out everything outside that range. \n
            Inputs: \n
                signal: a 1D numpy array to be filtered. \n
                highcut: the frequency (int or float) in Hz above which everything else will be filtered. \n
                lowcut: the frequency (int or float) in Hz below which everything else will be filtered. \n
                forder: the order (int) of the filter. \n
                fs: the sampling frequency (int or float) in Hz of the signal to be filtered. \n
            Outputs: \n
                filt_sig: a 1D numpy array representing the filtered signal. \n
            Sample call: \n
                filtered_signal = lowpass_filter(signal,lowcut,highcut,forder,fs)
        """
        
        win = signal.windows.tukey(len(self.data[:,0]),0.1)
        nyq = 0.5 * self.fs
        low = lowcut/nyq
        high = highcut/nyq
        
        b,a = signal.butter(forder, [low,high], btype='bandpass')
        
        # Apply the bandpass filter to all traces at once
        data_bandpass = signal.filtfilt(b, a, self.data * win[:, np.newaxis], axis=0)

        if inplace==True:
            self.data = data_bandpass
        else:
            return data_bandpass    

    def optimize_signature_window(self, trace_number, start_time_ms, end_time_ms,
                                  time_search_ms=2.0, length_search_ms=0.5,
                                  criterion='spectral_flatness', plot=False):
        """
        Automatically optimize the signature wavelet extraction window.
        
        This method searches around the specified time window to find the optimal
        wavelet extraction parameters based on various quality criteria.
        
        Parameters
        ----------
        trace_number : int
            Trace number to extract signature from (0-indexed)
        start_time_ms : float
            Initial start time in milliseconds (will be refined)
        end_time_ms : float
            Initial end time in milliseconds (will be refined)
        time_search_ms : float, optional
            Search range in milliseconds to shift the window. Default is 2.0 ms.
            The window will be tested at positions from -time_search_ms to +time_search_ms.
        length_search_ms : float, optional
            Search range for window length adjustment. Default is 0.5 ms.
            The window length will be tested from original ± length_search_ms.
        criterion : str, optional
            Optimization criterion. Default is 'spectral_flatness'.
            Options:
            - 'spectral_flatness': Maximize spectral flatness (whitest spectrum)
            - 'energy': Maximize RMS energy (strongest signal)
            - 'kurtosis': Maximize kurtosis (most impulsive/spiky)
            - 'bandwidth': Maximize bandwidth (broadest spectrum)
        plot : bool, optional
            If True, displays diagnostic plots. Default is False.
            
        Returns
        -------
        dict
            Dictionary with optimized parameters:
            - 'start_time_ms': Optimized start time
            - 'end_time_ms': Optimized end time
            - 'quality_score': Quality metric value
            - 'signature': Extracted signature wavelet
            
        Examples
        --------
        >>> seis = Seismic(data, fs=50000, dx=0.4)
        >>> result = seis.optimize_signature_window(6415, 30.5, 31.7)
        >>> print(f"Optimized window: {result['start_time_ms']:.2f} - {result['end_time_ms']:.2f} ms")
        >>> # Use optimized parameters for deconvolution
        >>> seis.signature_deconvolution(6415, result['start_time_ms'], result['end_time_ms'])
        """
        
        n_samples, n_traces = self.data.shape
        dt = 1.0 / self.fs
        
        # Validate trace number
        if trace_number < 0 or trace_number >= n_traces:
            raise ValueError(f"trace_number {trace_number} out of bounds: 0-{n_traces-1}")
        
        # Convert times to samples
        center_time_ms = (start_time_ms + end_time_ms) / 2.0
        initial_length_ms = end_time_ms - start_time_ms
        
        # Define search grid
        time_shifts = np.linspace(-time_search_ms, time_search_ms, 21)  # 21 positions
        length_adjustments = np.linspace(-length_search_ms, length_search_ms, 11)  # 11 lengths
        
        best_score = -np.inf
        best_params = None
        results = []
        
        print(f"Optimizing signature window for trace {trace_number}...")
        print(f"  Initial window: {start_time_ms:.2f} - {end_time_ms:.2f} ms ({initial_length_ms:.2f} ms long)")
        print(f"  Search: ±{time_search_ms:.2f} ms position, ±{length_search_ms:.2f} ms length")
        print(f"  Criterion: {criterion}")
        
        # Search over all combinations
        for time_shift in time_shifts:
            for length_adj in length_adjustments:
                # Calculate new window
                new_center = center_time_ms + time_shift
                new_length = initial_length_ms + length_adj
                new_start_ms = new_center - new_length / 2.0
                new_end_ms = new_center + new_length / 2.0
                
                # Convert to samples
                new_start_samp = int(new_start_ms * self.fs / 1000.0)
                new_end_samp = int(new_end_ms * self.fs / 1000.0)
                
                # Check bounds
                if new_start_samp < 0 or new_end_samp >= n_samples or new_start_samp >= new_end_samp:
                    continue
                
                # Extract signature
                sig = self.data[new_start_samp:new_end_samp, trace_number]
                
                if len(sig) < 5:  # Minimum length check
                    continue
                
                # Calculate quality metric
                score = self._calculate_wavelet_quality(sig, criterion)
                
                results.append({
                    'start_ms': new_start_ms,
                    'end_ms': new_end_ms,
                    'length_ms': new_length,
                    'score': score,
                    'signature': sig
                })
                
                if score > best_score:
                    best_score = score
                    best_params = {
                        'start_time_ms': new_start_ms,
                        'end_time_ms': new_end_ms,
                        'quality_score': score,
                        'signature': sig
                    }
        
        if best_params is None:
            raise ValueError("Could not find valid signature window in search range")
        
        print(f"  ✓ Optimized window: {best_params['start_time_ms']:.2f} - {best_params['end_time_ms']:.2f} ms")
        print(f"  Quality score: {best_score:.4f}")
        print(f"  Improvement: {best_params['start_time_ms'] - start_time_ms:+.2f} ms start, "
              f"{best_params['end_time_ms'] - end_time_ms:+.2f} ms end")
        
        # Optional plotting
        if plot:
            self._plot_signature_optimization(results, best_params, trace_number, 
                                             start_time_ms, end_time_ms, criterion)
        
        return best_params
    
    def _calculate_wavelet_quality(self, signature, criterion):
        """Calculate quality metric for a signature wavelet."""
        
        if criterion == 'energy':
            # RMS energy
            return np.sqrt(np.mean(signature**2))
        
        elif criterion == 'kurtosis':
            # Kurtosis (fourth moment) - measures "spikiness"
            if len(signature) < 4:
                return -np.inf
            mean = np.mean(signature)
            std = np.std(signature)
            if std == 0:
                return -np.inf
            return np.mean(((signature - mean) / std)**4)
        
        elif criterion == 'spectral_flatness':
            # Spectral flatness (ratio of geometric to arithmetic mean of power spectrum)
            # Higher = flatter/whiter spectrum = better for deconvolution
            sig_fft = np.fft.rfft(signature)
            power = np.abs(sig_fft)**2
            power = power[power > 0]  # Remove zeros for log
            if len(power) == 0:
                return -np.inf
            geometric_mean = np.exp(np.mean(np.log(power)))
            arithmetic_mean = np.mean(power)
            if arithmetic_mean == 0:
                return -np.inf
            return geometric_mean / arithmetic_mean
        
        elif criterion == 'bandwidth':
            # Effective bandwidth (based on -3dB points)
            sig_fft = np.fft.rfft(signature)
            power = np.abs(sig_fft)**2
            max_power = np.max(power)
            if max_power == 0:
                return -np.inf
            power_normalized = power / max_power
            # Find frequencies above -3dB (0.5 in linear scale)
            above_half_power = power_normalized > 0.5
            if not np.any(above_half_power):
                return -np.inf
            bandwidth = np.sum(above_half_power)  # Count of frequency bins
            return float(bandwidth)
        
        else:
            raise ValueError(f"Unknown criterion: {criterion}")
    
    def _plot_signature_optimization(self, results, best_params, trace_number, 
                                     original_start, original_end, criterion):
        """Plot diagnostic information about signature optimization."""
        import matplotlib.pyplot as plt
        
        # Convert results to arrays for plotting
        starts = np.array([r['start_ms'] for r in results])
        ends = np.array([r['end_ms'] for r in results])
        lengths = np.array([r['length_ms'] for r in results])
        scores = np.array([r['score'] for r in results])
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle(f'Signature Window Optimization - Trace {trace_number}', fontsize=14, fontweight='bold')
        
        # 1. Score vs start time
        ax = axes[0, 0]
        scatter = ax.scatter(starts, scores, c=lengths, cmap='viridis', s=30, alpha=0.6)
        ax.axvline(original_start, color='r', linestyle='--', alpha=0.5, label='Original')
        ax.axvline(best_params['start_time_ms'], color='g', linestyle='-', linewidth=2, label='Optimized')
        ax.set_xlabel('Start time (ms)')
        ax.set_ylabel(f'Quality score ({criterion})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax, label='Window length (ms)')
        
        # 2. Score vs window length
        ax = axes[0, 1]
        ax.scatter(lengths, scores, c=starts, cmap='plasma', s=30, alpha=0.6)
        original_length = original_end - original_start
        best_length = best_params['end_time_ms'] - best_params['start_time_ms']
        ax.axvline(original_length, color='r', linestyle='--', alpha=0.5, label='Original')
        ax.axvline(best_length, color='g', linestyle='-', linewidth=2, label='Optimized')
        ax.set_xlabel('Window length (ms)')
        ax.set_ylabel(f'Quality score ({criterion})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 3. 2D heatmap of scores
        ax = axes[0, 2]
        # Create grid for heatmap
        unique_starts = np.unique(starts)
        unique_lengths = np.unique(lengths)
        score_grid = np.full((len(unique_lengths), len(unique_starts)), np.nan)
        for r in results:
            i = np.where(unique_lengths == r['length_ms'])[0][0]
            j = np.where(unique_starts == r['start_ms'])[0][0]
            score_grid[i, j] = r['score']
        
        im = ax.imshow(score_grid, aspect='auto', origin='lower', cmap='RdYlGn',
                      extent=[unique_starts[0], unique_starts[-1], 
                             unique_lengths[0], unique_lengths[-1]])
        ax.plot(original_start, original_length, 'r*', markersize=15, label='Original')
        ax.plot(best_params['start_time_ms'], best_length, 'g*', markersize=20, label='Optimized')
        ax.set_xlabel('Start time (ms)')
        ax.set_ylabel('Window length (ms)')
        ax.set_title('Score Heatmap')
        ax.legend()
        plt.colorbar(im, ax=ax, label='Score')
        
        # 4. Original wavelet
        ax = axes[1, 0]
        dt = 1.0 / self.fs
        start_samp_orig = int(original_start * self.fs / 1000.0)
        end_samp_orig = int(original_end * self.fs / 1000.0)
        sig_orig = self.data[start_samp_orig:end_samp_orig, trace_number]
        time_orig = np.arange(len(sig_orig)) * dt * 1000 + original_start
        ax.plot(time_orig, sig_orig, 'r-', linewidth=1.5, label='Original window')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Amplitude')
        ax.set_title('Original Signature')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 5. Optimized wavelet
        ax = axes[1, 1]
        sig_opt = best_params['signature']
        time_opt = np.arange(len(sig_opt)) * dt * 1000 + best_params['start_time_ms']
        ax.plot(time_opt, sig_opt, 'g-', linewidth=1.5, label='Optimized window')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Amplitude')
        ax.set_title('Optimized Signature')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 6. Frequency comparison
        ax = axes[1, 2]
        fft_orig = np.fft.rfft(sig_orig)
        fft_opt = np.fft.rfft(sig_opt)
        freqs_orig = np.fft.rfftfreq(len(sig_orig), dt)
        freqs_opt = np.fft.rfftfreq(len(sig_opt), dt)
        ax.semilogy(freqs_orig/1000, np.abs(fft_orig)**2, 'r-', linewidth=1, alpha=0.7, label='Original')
        ax.semilogy(freqs_opt/1000, np.abs(fft_opt)**2, 'g-', linewidth=1.5, label='Optimized')
        ax.set_xlabel('Frequency (kHz)')
        ax.set_ylabel('Power Spectrum')
        ax.set_title('Frequency Content Comparison')
        ax.set_xlim([0, 20])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = f'signature_optimization_trace{trace_number}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved diagnostic plot: {filename}")
        plt.show()
    
    def signature_deconvolution(self, trace_number, start_time_ms, end_time_ms,
                                method='wiener', epsilon=0.01, prewhiten=True, 
                                prewhiten_percent=1.0, auto_optimize=False,
                                optimize_criterion='spectral_flatness', inplace=True):
        """
        Deconvolve seismic data using a signature wavelet extracted from a specified trace.
        This improves vertical resolution by removing the source/system signature.
        
        The deconvolution is applied in the frequency domain and is fully vectorized
        for efficient processing of multiple traces.
        
        Parameters
        ----------
        trace_number : int
            Trace number to extract the signature wavelet from (0-indexed)
        start_time_ms : float
            Start time in milliseconds for signature extraction window
        end_time_ms : float
            End time in milliseconds for signature extraction window
        method : str, optional
            Deconvolution method: 'wiener', 'spiking', or 'water-level'. 
            Default is 'wiener'.
            - 'wiener': Wiener deconvolution (frequency domain, stable)
            - 'spiking': Spiking deconvolution (assumes minimum phase)
            - 'water-level': Water-level deconvolution (stabilized inverse)
        epsilon : float, optional
            Stabilization parameter (noise level). Default is 0.01 (1%).
            Typical range: 0.001-0.1
        prewhiten : bool, optional
            Apply pre-whitening to signature. Default is True (recommended).
        prewhiten_percent : float, optional
            Pre-whitening percentage (0-100). Default is 1.0%.
            Typical range: 0.5-5.0
        auto_optimize : bool, optional
            If True, automatically optimizes the signature extraction window.
            Default is False. Uses optimize_signature_window() internally.
        optimize_criterion : str, optional
            Criterion for automatic optimization (if auto_optimize=True).
            Default is 'spectral_flatness'.
            Options: 'spectral_flatness', 'energy', 'kurtosis', 'bandwidth'
        inplace : bool, optional
            If True, modifies self.data. If False, returns deconvolved data.
            Default is True.
            
        Returns
        -------
        numpy.ndarray or None
            If inplace=False, returns deconvolved data array.
            If inplace=True, modifies self.data in place and returns None.
            
        Examples
        --------
        >>> seis = Seismic(data, fs=1000, dx=1.0)
        >>> # Manual window specification
        >>> seis.signature_deconvolution(50, 10, 30, method='wiener', epsilon=0.01)
        >>> 
        >>> # With automatic window optimization
        >>> seis.signature_deconvolution(50, 10, 30, auto_optimize=True)
        """
        
        # Auto-optimize window if requested
        if auto_optimize:
            print("Auto-optimizing signature window...")
            opt_result = self.optimize_signature_window(
                trace_number, start_time_ms, end_time_ms,
                criterion=optimize_criterion, plot=False
            )
            start_time_ms = opt_result['start_time_ms']
            end_time_ms = opt_result['end_time_ms']
            print(f"Using optimized window: {start_time_ms:.2f} - {end_time_ms:.2f} ms")
        
        # Validate inputs
        n_samples, n_traces = self.data.shape
        
        if trace_number < 0 or trace_number >= n_traces:
            raise ValueError(f"trace_number {trace_number} is out of bounds. Valid range: 0-{n_traces-1}")
        
        # Convert time to sample indices (fs is in Hz, so convert ms to seconds first)
        start_sample = int(start_time_ms * self.fs / 1000.0)
        end_sample = int(end_time_ms * self.fs / 1000.0)
        
        # Calculate max time in ms
        max_time_ms = (n_samples / self.fs) * 1000.0
        
        # Validate sample indices
        if start_sample < 0 or start_sample >= n_samples:
            raise ValueError(f"start_time_ms {start_time_ms} is out of bounds. Valid range: 0-{max_time_ms:.2f} ms")
        if end_sample <= start_sample or end_sample > n_samples:
            raise ValueError(f"end_time_ms {end_time_ms} is invalid. Must be > {start_time_ms} and <= {max_time_ms:.2f} ms")
        
        # Extract signature wavelet from specified trace and time window
        signature = self.data[start_sample:end_sample, trace_number].copy()
        n_sig = len(signature)
        
        if n_sig == 0:
            raise ValueError(f"Extracted signature has zero length. Check time window: {start_time_ms}-{end_time_ms} ms")
        
        n_data = self.data.shape[0]
        
        # Pad signature to match data length for frequency domain operations
        signature_padded = np.zeros(n_data)
        signature_padded[:n_sig] = signature
        
        # Apply pre-whitening to signature (recommended for stability)
        if prewhiten:
            whitening_factor = prewhiten_percent / 100.0
            max_sig_amp = np.max(np.abs(signature))
            if max_sig_amp > 0:
                signature_padded[:n_sig] = signature + whitening_factor * max_sig_amp
            else:
                print("Warning: Signature has zero amplitude, skipping pre-whitening")
        
        # Transform to frequency domain (vectorized for all traces)
        sig_fft = np.fft.fft(signature_padded)
        data_fft = np.fft.fft(self.data, axis=0)
        
        # Apply deconvolution based on method
        if method == 'wiener':
            # Wiener deconvolution: optimal in least-squares sense
            sig_power = np.abs(sig_fft)**2
            # Use constant stabilization (noise power estimate)
            # This provides uniform regularization across all frequencies
            noise_power = epsilon * np.mean(sig_power)  # Estimate of noise power
            wiener_filter = np.conj(sig_fft) / (sig_power + noise_power)
            # Apply to all traces (broadcasting)
            deconvolved_fft = data_fft * wiener_filter[:, np.newaxis]
            
        elif method == 'spiking':
            # Spiking deconvolution: assumes minimum phase
            sig_power = np.abs(sig_fft)**2
            # Inverse filter with water-level stabilization
            inverse_filter = 1.0 / (sig_fft + epsilon * np.max(np.abs(sig_fft)))
            deconvolved_fft = data_fft * inverse_filter[:, np.newaxis]
            
        elif method == 'water-level':
            # Water-level deconvolution: prevents division by small values
            water_level = epsilon * np.max(np.abs(sig_fft))
            # Create stabilized inverse
            sig_fft_stab = np.where(np.abs(sig_fft) < water_level, 
                                   water_level * np.exp(1j * np.angle(sig_fft)), 
                                   sig_fft)
            inverse_filter = 1.0 / sig_fft_stab
            deconvolved_fft = data_fft * inverse_filter[:, np.newaxis]
            
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'wiener', 'spiking', or 'water-level'.")
        
        # Transform back to time domain
        deconvolved_data = np.real(np.fft.ifft(deconvolved_fft, axis=0))
        
        # Update or return
        if inplace:
            self.data = deconvolved_data
            print(f"Signature deconvolution applied using {method} method")
            print(f"Signature extracted from trace {trace_number}, time {start_time_ms}-{end_time_ms} ms")
        else:
            return deconvolved_data

    def top_mute(self, horizon, shift_ms=0.0, taper_ms=0.0, inplace=True):
        """
        Apply a top mute to zero out amplitudes above a specified horizon.
        
        This is commonly used to remove noise or unwanted data above the seabed
        or first breaks. The mute can be shifted and optionally tapered to avoid
        sharp discontinuities.
        
        Parameters
        ----------
        horizon : float or array-like
            Mute time(s) in milliseconds.
            - float: constant mute time applied to all traces
            - array: time value for each trace (length must equal number of traces)
        shift_ms : float, optional
            Shift the mute boundary in milliseconds. Default is 0.0.
            - Positive values: mute MORE data (extend mute zone deeper)
            - Negative values: mute LESS data (preserve more signal)
            Example: shift_ms=-1.0 mutes everything above (horizon - 1 ms)
        taper_ms : float, optional
            Length of cosine taper in milliseconds at the mute boundary.
            Default is 0.0 (hard cutoff).
            The taper smoothly ramps from 0 (muted) to 1 (preserved) over this duration.
        inplace : bool, optional
            If True, modifies self.data in place. If False, returns muted data.
            Default is True.
            
        Returns
        -------
        numpy.ndarray or None
            If inplace=False, returns muted data array.
            If inplace=True, modifies self.data in place and returns None.
            
        Examples
        --------
        >>> seis = Seismic(data, fs=50000, dx=0.4)
        >>> # Mute above a constant time
        >>> seis.top_mute(horizon=10.0)
        >>> 
        >>> # Mute 1 ms above seabed pick (from Segy_edit object)
        >>> seis.top_mute(horizon=s.seabed_pick, shift_ms=-1.0)
        >>> 
        >>> # Mute with 2 ms cosine taper
        >>> seis.top_mute(horizon=s.seabed_pick, shift_ms=-1.0, taper_ms=2.0)
        
        Notes
        -----
        - Horizon times represent the mute boundary (data above is zeroed)
        - With shift_ms > 0: mute extends DEEPER (more conservative)
        - With shift_ms < 0: mute is SHALLOWER (preserves more data)
        - Taper is applied BELOW the mute boundary (in preserved zone)
        """
        
        n_samples, n_traces = self.data.shape
        dt = 1.0 / self.fs  # sampling interval in seconds
        
        # Handle horizon input
        if isinstance(horizon, (int, float)):
            # Constant horizon for all traces
            horizon_times = np.full(n_traces, horizon)
        else:
            # Array of horizon times
            horizon_times = np.asarray(horizon)
            if len(horizon_times) != n_traces:
                raise ValueError(f"Horizon length ({len(horizon_times)}) must match number of traces ({n_traces})")
        
        # Check for NaN values and handle them
        if np.any(np.isnan(horizon_times)):
            print(f"Warning: {np.sum(np.isnan(horizon_times))} NaN values in horizon. These traces will not be muted.")
        
        # Apply shift
        adjusted_horizon = horizon_times + shift_ms
        
        # Convert times from milliseconds to sample indices
        horizon_samples = (adjusted_horizon * 1e-3 / dt).astype(int)
        
        # Calculate taper length in samples
        taper_samples = int(taper_ms * 1e-3 / dt)
        
        # Create muted data (copy to avoid modifying original if inplace=False)
        if inplace:
            muted_data = self.data
        else:
            muted_data = self.data.copy()
        
        # Apply mute trace by trace
        for trace_idx in range(n_traces):
            horizon_samp = horizon_samples[trace_idx]
            
            # Skip if horizon is NaN
            if np.isnan(horizon_times[trace_idx]):
                continue
            
            # Skip if horizon is outside data range
            if horizon_samp < 0:
                # Horizon is before data start - mute entire trace
                muted_data[:, trace_idx] = 0.0
                continue
            elif horizon_samp >= n_samples:
                # Horizon is after data end - don't mute anything
                continue
            
            # Zero out data above horizon
            muted_data[:horizon_samp, trace_idx] = 0.0
            
            # Apply taper if requested
            if taper_ms > 0 and taper_samples > 0:
                # Define taper zone: from horizon_samp to horizon_samp + taper_samples
                taper_end = min(horizon_samp + taper_samples, n_samples)
                taper_length = taper_end - horizon_samp
                
                if taper_length > 0:
                    # Create cosine taper (0 at top, 1 at bottom)
                    # Using (1 - cos) / 2 for smooth 0 to 1 transition
                    taper = (1.0 - np.cos(np.linspace(0, np.pi, taper_length))) / 2.0
                    muted_data[horizon_samp:taper_end, trace_idx] *= taper
        
        # Print summary
        valid_horizons = ~np.isnan(horizon_times)
        if np.any(valid_horizons):
            min_horizon = np.nanmin(adjusted_horizon)
            max_horizon = np.nanmax(adjusted_horizon)
            mean_horizon = np.nanmean(adjusted_horizon)
            print(f"Top mute applied:")
            print(f"  Horizon range: {min_horizon:.2f} - {max_horizon:.2f} ms (mean: {mean_horizon:.2f} ms)")
            if shift_ms != 0:
                print(f"  Shift: {shift_ms:+.2f} ms")
            if taper_ms > 0:
                print(f"  Taper: {taper_ms:.2f} ms ({taper_samples} samples)")
        
        if not inplace:
            return muted_data

    def fk_spectrum(self, pad_t=2, pad_x=2):
        '''
        Calculates the 2D Fourier transform, for seismic data this results in
        the data being transformed into the F-K domain. Zero-padding of the input 
        data is performed for both axis before the transform is applied.
        
        Parameters
        ----------
        pad_t : int, optional
            Time axis zero-padding factor. The default is 2.
        pad_x : in, optional
            Spatial (trace) axis zero-adding factor. The default is 2.

        Returns
        -------
        data_fk : 2D Array of floats
            2D fourier transformed data (F-K domain)
        f : 1D Array of floats
            Frequency axis vecor belonging to data_fk
        kx : 1D Array of floats
            Wavenumber axis vecor belonging to data_fk
        '''
        
        nt,nx = np.shape(self.data)
        
        nt_fft = pad_t*nt #padding for time axis
        nx_fft = pad_x*nx #padding for trace axis
       
        
        fk_full = fft2(self.data, s=(nt_fft, nx_fft))/len(self.data)  
        
        fk_shift = fftshift((np.abs(fk_full)))
        data_fk = fk_shift[int(len(fk_full)/2):,:]
        
        f =  np.linspace(-0.5,0.5, num=2*fk_full.shape[0])*self.fs
        kx = np.linspace(-0.5,0.5, num=2*fk_full.shape[1])/self.dx
        
        self.data_fk = data_fk
        self.f = f
        self.kx = kx

    
    def plot_fk(self, clip=1, kwin=0, fwin=0, pad_t=2, pad_x=2, outfile=None):  #NEEDS PROPER QC
        '''
        Plot of FK-spectrum and optionally save as figure

        Parameters
        ----------
        gain : float, optional
            Gains up the amplitudes in the F-K spectrum to
            tweak optimal visualisation. The default is 1.
        fwin : float, optional
            frequency window: zoom from 0 to this amount of Hz .
            The default is 0 and plots all available frequencies
        pad_t : int, optional
            Time axis zero-padding factor. The default is 2.
        pad_x : in, optional
            Spatial (trace) axis zero-adding factor. The default is 2.
        outfile : str, optional
            To save the figure, set outfile to e.g.
            outfile=r'C:/myfolder/fk-spectrum.png 
            The default is None.

        Returns
        -------
        ax: image.AxesImage
            this can be used for customized plotting e.g. in subplots

        '''
        
        try:
            vmin = clip*np.min(self.data_fk)
            vmax = clip*np.max(self.data_fk)
        except AttributeError:
            print("'Seismic' object has no attribute 'data_fk', calculate fk_spectrum first")
            sys.exit(0)
        
        
        lognorm = colors.LogNorm(vmin=vmin, vmax=vmax)
        
        if fwin==0:
            extent = [self.kx[0], self.kx[-1], np.median(self.f), self.f[-1]]
            data_fk_plot = self.data_fk[0: int(self.data_fk.shape[0]),:]
            
        else: #very hacky, improve
            idx = np.argmax(self.f >= fwin)
            idx0 = np.argmax(self.f >= 0)
            idx_pos = idx - idx0
            idx_data = idx_pos * (len(self.data_fk)/(0.5*len(self.f)))
            
            extent = [self.kx[0], self.kx[-1], np.median(self.f), self.f[idx]]
            
            vmin = clip*np.min(self.data_fk[0: int(idx_data),:])
            vmax = clip*np.max(self.data_fk[0: int(idx_data),:])
            lognorm = colors.LogNorm(vmin=vmin, vmax=vmax)
            
            data_fk_plot = self.data_fk[0: int(idx_data),:]
            
        if kwin != 0:
            idx_min = np.argmax(self.kx >= -kwin) -1
            idx_max = np.argmax(self.kx > kwin)
            
            extent = [self.kx[idx_min], self.kx[idx_max], np.median(self.f), self.f[idx]]
            data_fk_plot = data_fk_plot[:, int(idx_min/2):int(idx_max/2)]
            
            
       
        cmap = 'jet'
        fig, ax = plt.subplots(1,1, figsize=(4, 4))
        
        
        im =  ax.imshow(data_fk_plot, cmap=cmap, extent=extent, aspect='auto',  norm=lognorm, interpolation='kaiser', origin='lower')
        ax.set_title('F-K spectrum, Fs = {} Hz, dx = {} m'.format(self.fs, np.round(self.dx,2)))
        ax.set_ylabel('Frequency (Hz)')
        ax.set_xlabel('Wave number (1/m)')
        ax.invert_xaxis()
            
        cbar = fig.colorbar(im)
        cbar.set_label("Amplitude (lognorm)", rotation=270,labelpad=20)
        
        plt.tight_layout()
        
        try:
            plt.savefig(outfile, dpi=300)
        except AttributeError:
            print("Set 'outpath' to save the frequency-wavenumber figure")
            pass
        
        return ax
        
            
    def plot_time_series(self, data=None, clip=1, win=0, outfile=None):
        
        '''
        Plot of timeseries and optionally save as figure

        Parameters
        ----------
        clip : float, optional
            Percentage clip of the min and max value used for
            visualisation. The default is 1.
        win : float, optional
            time window: zoom from start of record to this amount of seconds 
   
        outfile : str, optional
            To save the figure, set outfile to e.g.
            outfile=r'C:/myfolder/time_series.png 
            The default is None.

        Returns
        -------
        ax: image.AxesImage
            this can be used for customized plotting e.g. in subplots
        
        '''
        
        if data is None:
            data = self.data
           
       
              
        fig, ax = plt.subplots(1,figsize=(4, 4))
        fig.suptitle('Time series, Fs = {} Hz, dx = {} m'.format(self.fs, np.round(self.dx,2)))
        
        if win==0:
            extent = [0, data.shape[1]*self.dx, data.shape[0]/self.fs, 0]
            
            amin = clip*data.min()
            amax = clip*data.max()
                    
            im = ax.imshow(data, cmap='Greys', aspect='auto', extent=extent, vmin=amin, vmax=amax)
        
        else:
           extent = [0, data.shape[1]*self.dx, win, 0]
           samples = int(win*self.fs)
           
           amin = clip*data[:samples,:].min()
           amax = clip*data[:samples,:].max()
                 
           im = ax.imshow(data[:samples,:], cmap='Greys', aspect='auto', interpolation='kaiser', extent=extent, vmin=amin, vmax=amax)
            
        
        ax.set_ylabel('Time (seconds)')
        ax.set_xlabel('Distance (m)')
        
        cbar = fig.colorbar(im)
        cbar.set_label("Amplitude", rotation=270,labelpad=20)
        cbar.formatter.set_powerlimits((0, 0))
        plt.tight_layout()
                
        try:
            plt.savefig(outfile, dpi=300)
        except AttributeError:
            print("Set 'outpath' to save the timeseries figure")
            pass
        
        return ax
        
    def difference_plot(self, data1, data2, clip=1, normalize=False, win=0, invert=False, outfile=None):
        
        '''
        Plot the difference between two datasets to show the effect
        of a specific processing step
        
                
        Parameters
        ----------
        data1: 2d array
            First dataset (e.g. seismic shot)
        data2: 2d array
            Second dataset (e.g. filtered seismic shot)
            that is to be extracted from the first
        clip : float, optional
            Percentage clip of the min and max value used for
            visualisation. The default is 1.
        normalize: bool, optional
            set the min and max value for the plot based on the 
            resulting difference plot (True) or to the first dataset (False).
            The default is False
        win : float, optional
            time window: zoom from start of record to this amount of seconds 
        invert: bool, optional
            invert the polarit of the resulting difference plot. 
            Default is False.
        outfile : str, optional
            To save the figure, set outfile to e.g.
            outfile=r'C:/myfolder/time_series.png 
            The default is None.

        Returns
        -------
        ax: image.AxesImage
            this can be used for customized plotting e.g. in subplots
        
        '''
          
        difference = data1 - data2
        
        if invert==True:
            difference = -difference
        
        if normalize==True:
            dmin = clip*difference.min()
            dmax = clip*difference.max()
            
        else: 
           dmin = clip*data1.min()
           dmax = clip*data1.max()
              
        fig, ax = plt.subplots(1,figsize=(4, 4))
        fig.suptitle('Time series, Fs = {} Hz, dx = {} m'.format(self.fs, np.round(self.dx,2)))
        
        if win==0:
            extent = [0, data1.shape[1]*self.dx, data1.shape[0]/self.fs, 0]
            im = ax.imshow(difference, cmap='Greys', aspect='auto', extent=extent, vmin=dmin, vmax=dmax)
        
        else:
           extent = [0, data1.shape[1]*self.dx, win, 0]
           samples = int(win*self.fs)
           im = ax.imshow(difference[:samples,:], cmap='Greys', aspect='auto', extent=extent, vmin=dmin, vmax=dmax)
            
        
        ax.set_ylabel('Time (seconds)')
        ax.set_xlabel('Distance (m)')
        
        cbar = fig.colorbar(im)
        cbar.set_label("Amplitude difference", rotation=270,labelpad=20)
        plt.tight_layout()
                
        try:
            plt.savefig(outfile, dpi=300)
        except AttributeError:
            print("Set 'outpath' to save the timeseries difference figure")
            
        return ax
        
    def plot_histogram(self, bins=256, logscale=False, outfile=None):
        
        '''
        Plot of histogram of the amplitudes of the data

        Parameters
        ----------
        bins : int, optional
           Sets the amount of bins to use for the
           histogram. The default is 256.
        logscale : bool, optional
            Sets the y-axis of the histogram as
            logscale. The default is False

        Returns
        -------
        hist: tuple with
            - n: values of histogram bins
            - bins: the edges of the bins
            - patches: .BarContainer object  
        
        This is output for customized plotting
                    
        '''
            
        plt.figure(figsize=(12,4))   
        hist = plt.hist(self.data.ravel(), bins=bins, range=[self.data.min().astype('float32'), self.data.max().astype('float32')])
        plt.title("Histogram of Amplitudes ({} bins)".format(bins))
        plt.xlabel('Signal amplitude ({})'.format(self.byte_format))
        plt.ylabel('Frequency of occurence')
        
        plt.grid()
        
        if logscale==True:
            plt.yscale('log')
            plt.ylabel('Frequency of occurence (log)')
            
        try:
            plt.savefig(outfile, dpi=300)
        except AttributeError:
            print("Set 'outpath' to save the histogram figure")
        
        return hist
                
    def plot_psd(self, clip=1, fwin=0, outfile=None):
        
        '''
        Plot of power spectral density and optionally save as figure

        Parameters
        ----------
        clip : float, optional
            Percentage clip of the min and max value used for
            visualisation. The default is 1.
        fwin : float, optional
           frequency window: zoom from 0 to this amount of Hz .
           The default is 0 and plots all available frequencies
        outfile : str, optional
            To save the figure, set outfile to e.g.
            outfile=r'C:/myfolder/psd_plot.png 
            The default is None.

        Returns
        -------
        ax: image.AxesImage
            this can be used for customized plotting e.g. in subplots
        
        '''
            
        self.f_psd, self.Pxx_den = sgn.welch(self.data.transpose() ,self.fs, nperseg=None, noverlap=None, nfft=self.data.shape[0])
              
        cmap='jet'
        
        fig, ax = plt.subplots(1, figsize=(4, 4))
        fig.suptitle('Power spectral density, Fs = {} Hz, dx = {} m'.format(self.fs, np.round(self.dx,2)))
       
        if fwin==0:
            extent = [0, self.data.shape[1]*self.dx, 0, self.f_psd.max()]
            
            pmin = clip*self.Pxx_den.min()
            pmax = clip*self.Pxx_den.max()
            
            lognorm = colors.LogNorm(vmin=pmin, vmax=pmax)
            im = ax.imshow(self.Pxx_den.transpose(), cmap=cmap, extent=extent, aspect='auto',  norm=lognorm, interpolation='kaiser', origin='lower')
            
        else:
            fidx = np.argwhere(self.f_psd <= fwin).max()
            extent = [0, self.data.shape[1]*self.dx, 0, self.f_psd[fidx]]
            
            pmin = clip*self.Pxx_den[:, :fidx].min()
            pmax = clip*self.Pxx_den[:, :fidx].max()
                        
            lognorm = colors.LogNorm(vmin=pmin, vmax=pmax)
            im = ax.imshow(self.Pxx_den[:, :fidx].transpose(), cmap=cmap, extent=extent, aspect='auto',  norm=lognorm, interpolation='kaiser', origin='lower')
        
        ax.set_ylabel('Frequency (Hz)')
        ax.set_xlabel('Distance')
        
        
        cbar = fig.colorbar(im)
        cbar.set_label("Amplitude V**2/Hz (lognorm)", rotation=270,labelpad=20)
        plt.tight_layout()
        
        
        try:
            plt.savefig(outfile, dpi=300)
        except AttributeError:
            print("Set 'outpath' to save the power spectral density figure")
            pass
        
        return ax
       

    
    
    
    
    
    
    
    
    
    