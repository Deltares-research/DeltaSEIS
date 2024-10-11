# -*- coding: utf-8 -*-
"""
Testscript for the deltaseis module

@author: nieboer
"""
#%%
from deltaseis_module import Seismic
import numpy as np
from obspy import read
import tdms_reader
import matplotlib.pyplot as plt
from pathlib import Path

#%%
def imshow_kwargs(ax):
    
    'Retrieve all the imshow keyword arguments from existing figure'
    
    im = ax.images[0]
    imshow_kwargs = {"X" : im.get_array(), "extent" : im.get_extent(), "aspect" : ax.get_aspect(),
                     "cmap" : im.get_cmap(), "interpolation" : im.get_interpolation(), 
                     "norm" : im.norm, "origin" : im.origin}
   
    return imshow_kwargs 

def time_psd_fk(time0, psd0, fk0, outfile):
    fig, axs = plt.subplots(1,3, figsize=(18,7))
    
    time_kwargs = imshow_kwargs(time0)
    psd_kwargs = imshow_kwargs(psd0)
    fk_kwargs = imshow_kwargs(fk0)
                          
    im0 = axs[0].imshow(**time_kwargs)
    axs[0].set_ylabel('Time (seconds)')
    axs[0].set_xlabel('Distance (m)')
    cbar = fig.colorbar(im0)
    cbar.set_label("Amplitude", rotation=270,labelpad=20)
    cbar.formatter.set_powerlimits((0, 0))
     
    im1 = axs[1].imshow(**psd_kwargs)
    axs[1].set_ylabel('Frequency (Hz)')
    axs[1].set_xlabel('Distance')
    cbar = fig.colorbar(im1)
    cbar.set_label("Amplitude V**2/Hz (lognorm)", rotation=270,labelpad=20)
    
    
    im2 = axs[2].imshow(**fk_kwargs)
    axs[2].set_ylabel('Frequency (Hz)')
    axs[2].set_xlabel('Wave number (1/m)')
    axs[2].invert_xaxis()
        
    cbar = fig.colorbar(im2)
    cbar.set_label("Amplitude (lognorm)", rotation=270,labelpad=20)
    plt.tight_layout(h_pad=10)
    plt.savefig(outfile)

#%%
#load hydrophone file
# sg2_file = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET43_CH5\30.dat")    #near offset
sg2_file = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET87_CH23\51.dat") #center offset
# sg2_file = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET148_CH48\39.dat")  #far offset
st = read(sg2_file, format="SEG2")

shot =[st[tr].data for tr in range(len(st))]  
data = np.vstack(shot).transpose()

props = st.traces[0].stats

fs = props.sampling_rate #Hz sampling frequency
dx = 5    #m, trace spacing

outpath = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\figures")
outfile = outpath/sg2_file.with_suffix('.png').name

#%%load data from iDAS tdsm of the helically wound cable


# tdms_path = Path("D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET43_CH5\Airgun__12m_UTC+0100_DST1_20230419_153240.453.tdms")     #near offset
# tdms_path = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET87_CH23\Airgun_UTC+0100_DST1_20230419_180448.574.tdms")          #center offset
tdms_path = Path(r"D:\Projects\Sweden\TEST_SWEDEN\data_selection\OFFSET148_CH48\Airgun_UTC+0100_DST1_20230419_170711.447.tdms")         #far offset

# t_start = 0.270  #ch05 shot
# t_start = 3.715    #ch23 shot
t_start = 2.450    #ch48 shot


helical_wind_factor = 1.15

tdms = tdms_reader.TdmsReader(tdms_path) 
data_tdms = tdms.get_data()

props_tdms = tdms.get_properties() 
fs_tdms = props_tdms.get("SamplingFrequency[Hz]")
dx_tdms = props_tdms.get('SpatialResolution[m]') * props_tdms.get('Fibre Length Multiplier') / helical_wind_factor


outfile_tdms = outpath/tdms_path.with_suffix('.png').name

 
#%%initiate instance with raw seismic data
seis = Seismic(data, fs, dx)

#calculate and plot amplitude histogram
# hist = seis.plot_histogram(logscale=True, bins=256)

#plot the time series of the shot
# time = seis.plot_time_series(clip=1, win=1.5)
# psd  = seis.plot_psd(clip=10000, fwin=50)

#calculate and plot the f-k spectrum
seis.fk_spectrum()
# fk = seis.plot_fk(fwin=50)

#plot results side by side
# time_psd_fk(time, psd, fk, outfile)


#%%optionally use a bandpass filter
seis.bandpass_filter(80, 300)
seis.agc_gain(seis.data_bandpass,agc_gate=1)

#%%apply agc gain to improve fk-spectrum
# seis.agc_gain(agc_gate=1)
seis_agc = Seismic(seis.data_agc, fs, dx)

#plot the time series of the shot
time_agc = seis_agc.plot_time_series(clip=0.1, win=0.3)
psd_agc  = seis_agc.plot_psd(clip=5, fwin=500)

#calculate and plot the f-k spectrum
seis_agc.fk_spectrum()
fk_agc = seis_agc.plot_fk(clip=10, fwin=500)


outfile_agc = Path(str(outfile).replace('.png', '_agc.png'))
#plot results side by side
time_psd_fk(time_agc, psd_agc, fk_agc, outfile_agc)


#%%initiate instance with raw seismic data

t_end = t_start + 1.5

ch_start = 120
ch_end = ch_start + 193 #max channel of helical cable 167*1.15

data_tdms_shot = data_tdms[int(t_start*fs_tdms):int(t_end*fs_tdms), ch_start:ch_end]

# seis_tdms = Seismic(data_tdms, fs_tdms, dx_tdms)     #to show the whole record
seis_tdms = Seismic(data_tdms_shot, fs_tdms, dx_tdms)  #to show the shot only, with 15 ms before first arrival

#calculate and plot amplitude histogram
# hist = seis_tdms.plot_histogram(logscale=True, bins=256)

#plot the time series of the shot
# time_tdms = seis_tdms.plot_time_series(clip=0.01)
# psd_tdms  = seis_tdms.plot_psd(clip=1, fwin=50)

#calculate and plot the f-k spectrum
seis_tdms.fk_spectrum()
# fk_tdms = seis_tdms.plot_fk(fwin=50)

#plot results side by side
# time_psd_fk(time_tdms, psd_tdms, fk_tdms, outfile_tdms)


#%%optionally apply bandpass filter
seis_tdms.bandpass_filter(40, 200, forder=1)
seis_tdms.agc_gain(seis_tdms.data_bandpass, agc_gate=1)

#%%apply agc gain to improve fk-spectrum
# seis_tdms.agc_gain(agc_gate=1)
seis_tdms_agc = Seismic(seis_tdms.data_agc, fs_tdms, dx_tdms) #load in banpass filtered and agc corrected data

#plot the time series of the shot
time_tdms_agc = seis_tdms_agc.plot_time_series(clip=0.5, win=0.3)


psd_tdms_agc  = seis_tdms_agc.plot_psd(fwin=500, clip=1)

#calculate and plot the f-k spectrum
seis_tdms_agc.fk_spectrum()
fk_tdms_agc = seis_tdms_agc.plot_fk(fwin=500, kwin=0.10, clip=10)


outfile_tdms_agc = Path(str(outfile_tdms).replace('.png', '_agc.png'))
#plot results side by side
time_psd_fk(time_tdms_agc, psd_tdms_agc, fk_tdms_agc, outfile_tdms_agc)

#%% coupling

seis_tdms.time_squared_gain()
seis_tdms_t2 = Seismic(seis_tdms.data_t2_gain, fs_tdms, dx_tdms)

seis_tdms.plot_time_series(clip=0.5, win=0.8)
seis_tdms_t2.plot_time_series(clip=0.1, win=0.8)


seis_tdms_t2.plot_psd(fwin=300)

uncoupled_distance = 82.3 #m
coupled_distance = 72.3 #m

uncoupled_channel = int(round(uncoupled_distance / seis_tdms.dx, 0))
coupled_channel = int(round(coupled_distance / seis_tdms.dx, 0))

print(seis_tdms_t2.Pxx_den.shape)
print(seis_tdms_t2.data.shape)

psd_uncoupled = seis_tdms_t2.Pxx_den[uncoupled_channel, :]
psd_coupled   = seis_tdms_t2.Pxx_den[coupled_channel, :]



plt.figure()
plt.plot(psd_uncoupled)
plt.plot(psd_coupled)

#%% different approach

uchan_time = seis_tdms_t2.data[:, uncoupled_channel]
cchan_time = seis_tdms_t2.data[:, coupled_channel]

#pad to total the next power 2 after twice (choice) the signal length
# N = 2**int(np.ceil(np.log2(2*len(uchan_time)))) 
N = int(len(uchan_time))


uchan_time = np.pad(uchan_time, N - len(uchan_time), 'constant')
cchan_time = np.pad(cchan_time, N - len(cchan_time), 'constant')

plt.figure(11)
plt.title('Time series')
plt.plot(uchan_time, label='Bad coupling')
plt.plot(cchan_time, label='Good coupling')
plt.legend()
plt.grid()

amp_uchan = np.abs(np.fft.fft(uchan_time))
amp_cchan = np.abs(np.fft.fft(cchan_time))



frequency_axis = np.fft.fftfreq(len(uchan_time), d=1/seis_tdms_t2.fs)

freq_win = 200
idx_fwin = np.argmin(np.abs(frequency_axis - freq_win))

plt.figure(22)
idx = int(len(frequency_axis)/2)

plt.title("Amplitude spectra")
plt.plot(frequency_axis[:idx][:idx_fwin], amp_uchan[:idx][:idx_fwin], label="Bad coupling")
plt.plot(frequency_axis[:idx][:idx_fwin], amp_cchan[:idx][:idx_fwin], label="Good coupling")
# plt.plot(frequency_axis, amp_cchan)
plt.legend()
plt.grid()





