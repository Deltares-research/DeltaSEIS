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
            self.data_bandpass = data_bandpass[:, 1:]
        else:
            return data_bandpass[:, 1:]    


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
       

    
    
    
    
    
    
    
    
    
    