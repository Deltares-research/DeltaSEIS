import numpy
import segpy
import math
from obspy.imaging.cm import viridis_white_r
import matplotlib.pylab as plt

#AFTER CHECK IF MASW AND MASW2 GIVE THE SAME RESULT (AND IF MASW IS FASTER), MOVE ALL THIS INTO BASE_SEISMIC CLASS (ILNCUDE PLOT)

def masw(Signal_masw,dx,freq_min,freq_max,vph_min,vph_max,fs,source_pos,y_offset):
    
    Data = Signal_masw     
    freq = np.arange(freq_min, freq_max, 0.2)
    vph = np.arange(vph_min, vph_max, 1)  
    dt = 1/fs     
    fmax = (1/dt)/2
    res = math.log2((2*(fmax/(freq[2]-freq[1]))))
    nfft = int(math.pow(2,res))    
    df = fmax/(nfft/2)  
    f = np.arange(0,fmax,df)
    st_fft = np.fft.fft(np.transpose(Data),nfft)
    st_fft = np.transpose(st_fft)
    Ry = st_fft[range(0,int(nfft/2)),:]
    
    offset = np.sqrt((((np.arange(0,len(Data[0,:]))*dx)+source_pos)**2 + (y_offset**2)))
    
    # finding frequency indixes
    
    freq_index = np.zeros((len(freq),1))
    
    for i in range(0,len(freq)):
        
        freq_index[i] = (np.abs(freq[i]-f)).argmin()
    
    PhaseInfo_r = np.angle(Ry[freq_index.astype(int),:])
    PhaseInfo = np.reshape(PhaseInfo_r,(len(freq),len(Data[0,:])))
    
    # Computing dispersion Image
    
    w = 2*np.pi*freq
    Energy = np.zeros((len(freq),len(vph)))
    
    PhaseVelSpectrum = np.zeros((len(freq),len(vph)))
    
    for ii in range (0,len(freq)):
        
        for jj in range(0,len(vph)):
            
            if (vph[jj]/freq[ii])<(dx/10):
                
                Energy[ii,jj] = 0
                
            else:
                Energy[ii,jj] += np.abs(np.sum((np.exp(1j*PhaseInfo[ii,:]))*np.exp(1j*(w[ii]*offset/vph[jj]))))
        PhaseVelSpectrum[ii,:] = np.abs(Energy[ii,:])     

    return(PhaseVelSpectrum,PhaseInfo) 

#%%

def masw2(data, dx, freq_min, freq_max, vph_min, vph_max, fs, source_pos, y_offset):

    ##AFTER QC MOVE THIS TO BASE SEISMCIC TO BE USED IN CLASS
    
    # Create frequency and phase velocity arrays
    freq = np.arange(freq_min, freq_max, 0.2)
    vph = np.arange(vph_min, vph_max, 1)
    
    # Compute constants
    dt = 1 / fs
    fmax = (1 / dt) / 2
    res = math.log2(2 * (fmax / (freq[2] - freq[1])))
    nfft = int(2 ** res)
    df = fmax / (nfft / 2)
    f = np.arange(0, fmax, df)
    
    # FFT of the input data
    st_fft = np.fft.fft(data.T, nfft).T
    Ry = st_fft[:nfft//2, :]
    
    # Compute offset array
    offset = np.sqrt((((np.arange(len(data[0])) * dx) + source_pos) ** 2) + (y_offset ** 2))
    
    # Find frequency indices
    freq_index = np.array([(np.abs(freq[i] - f)).argmin() for i in range(len(freq))])
    
    # Compute phase information
    PhaseInfo_r = np.angle(Ry[freq_index.astype(int), :])
    PhaseInfo = PhaseInfo_r.reshape(len(freq), len(data[0]))
    
    # Initialize Energy and PhaseVelSpectrum arrays
    w = 2 * np.pi * freq
    Energy = np.zeros((len(freq), len(vph)))
    
    # Use vectorized computation to avoid nested loops
    condition = (vph[:, np.newaxis] / freq) < (dx / 10)
    Energy = np.where(condition.T, 0, Energy)  # Apply condition to zero out energy
    
    # Vectorized calculation of Energy
    exp_phase = np.exp(1j * PhaseInfo)[:, np.newaxis, :]  # Shape (len(freq), 1, len(Data[0]))
    exp_offset = np.exp(1j * (w[:, np.newaxis] * offset / vph))  # Shape (len(freq), len(vph), len(Data[0]))
    
    # Summing across the third axis (the trace dimension) to update Energy
    Energy += np.abs(np.sum(exp_phase * exp_offset, axis=2))
    
    # Calculate PhaseVelSpectrum
    PhaseVelSpectrum = np.abs(Energy)
    
    return PhaseVelSpectrum, PhaseInfo
#%%

def plot_masw(data, phase_velocity_spectrum, outfile=None):

    fig, ax = plt.subplots(figsize=(15,5))
    plt.subplot(121)
    segypy.wiggle(sg_z,t=time,x=offset)
    plt.xlabel("Receiver coordinates [m]",fontsize=14)
    #plt.imshow(np.abs(sg_res),interpolation='kaiser', aspect='auto',cmap=viridis_white_r)
    plt.ylabel("Time [s]",fontsize=14)
    plt.title("Shotgather",fontsize=14)
    #plt.xlim(0,100)
    
    PhaseVelSpectrum = np.transpose(PhaseVelSpectrum)

    #PhaseVelSpectrum = PhaseVelSpectrum/np.max(PhaseVelSpectrum,axis=0)
        
    plt.subplot(122)
    plt.imshow(PhaseVelSpectrum,interpolation='kaiser', aspect='auto',cmap='jet',extent=[freq_min, freq_max, vph_max, vph_min])
    plt.contour(PhaseVelSpectrum,interpolation='kaiser', aspect='auto',cmap='jet',extent=[freq_min, freq_max, vph_min, vph_max])
    plt.xlabel("Frequency [Hz]",fontsize=14)
    plt.ylabel("Phase velocity [m/s]",fontsize=14)
    plt.title("Phase velocity Spectrum",fontsize=14)
    plt.xlim(freq_min,freq_max)
    plt.ylim(vph_min,vph_max)

    #%%
    freq = np.arange(freq_min, freq_max, 0.2)
    vph = np.arange(vph_min, vph_max, 1) 
    v_max =  np.array([vph[np.where(np.max(PhaseVelSpectrum,axis=0)[ik] == PhaseVelSpectrum[:,ik])] for ik in range(0,len(freq))])[:,0]

    plt.imshow(PhaseVelSpectrum,interpolation='kaiser', aspect='auto',cmap='jet',origin='lower',extent=[freq_min, freq_max, vph_min, vph_max])
    plt.contour(PhaseVelSpectrum,interpolation='kaiser', aspect='auto',cmap='jet',extent=[freq_min, freq_max, vph_min, vph_max])
    plt.xlabel("Frequency [Hz]",fontsize=14)
    plt.ylabel("Phase Velocity [m/sec]",fontsize=14)
    plt.title("Phase Velocity Spectrum",fontsize=14)
    plt.xlim(freq_min,freq_max)
    plt.ylim(vph_min,vph_max)

    data_c = np.loadtxt("theo_curve.txt")

    plt.plot(data_c[:,0],data_c[:,1])
    plt.plot(freq,v_max,'--',color='green')

    plt.show()

    if outfile is not None:
        plt.savefig(outfile)

    plt.close()
