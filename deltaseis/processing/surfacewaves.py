def masw(Signal_masw,dx,freq_min,freq_max,vph_min,vph_max,fs,source_pos,y_offset):

    
    import numpy as np
    import math
    
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