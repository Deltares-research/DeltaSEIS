# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 19:13:22 2023
 index out of bounds error at 172
 fix with edwin
 
 add docstrings

 do we stil need to export to sg2/dat as well? where do we use that (not vista)
 
 
@author: nieboer
"""
#%%
def header_geometry(ntraces, dx, shot_number, spi):

    import numpy as np
        
    geometry_dict = {}
    
    geometry_dict['zeros']          = [0]*ntraces
    geometry_dict['spn']            = [shot_number]*ntraces
    geometry_dict['ffid']           = [spn + 1000 for spn in geometry_dict['spn']]
    geometry_dict['source_number']  = [shot_number]*ntraces
    geometry_dict['channels']       = [channel + 1 for channel in range(ntraces)] 
    geometry_dict['x_rec']          = [spi + (i-1)*dx for i in geometry_dict['channels']]
    geometry_dict['y_rec']          = geometry_dict['zeros']
    geometry_dict['x_src']          = [(shot_number - 1)*spi]*ntraces
    geometry_dict['y_src']          = geometry_dict['zeros']
    geometry_dict['offset']         = [np.sqrt((y_rec - y_src)**2 + (x_rec - x_src)**2) for y_rec, y_src, x_rec, x_src in zip(geometry_dict['y_rec'],  geometry_dict['y_src'],  geometry_dict['x_rec'],  geometry_dict['x_src'])]
    geometry_dict['x_cdp']          = [0.5*(x_rec - x_src) + x_src for x_rec, x_src in zip(geometry_dict['x_rec'],geometry_dict['x_src'])]
    geometry_dict['y_cdp']          = geometry_dict['zeros']
    geometry_dict['cdp']            = [int(x_cdp/(0.5*dx)) for x_cdp in geometry_dict['x_cdp']]
    
    c0 = int( geometry_dict['x_rec'][0]/dx) #trick to let the first CDP number start at 1
    
    geometry_dict['cdp']        = [cdp - c0 + 1 for cdp in geometry_dict['cdp']]  #set first cdp number to 1   
    
    return geometry_dict

def export_sgy(data, fs, dx, shot_number, spi, out_sgy):
    import numpy as np
    import segyio
    
    data = np.float32(data)
    #data = np.ascontiguousarray(data, dtype=np.float32)
    print("\nData was converted to 32-bit floating point numbers\n")
    
    nsamples, ntraces = np.shape(data)
    recording_delay = 0
    scaler = -100
    
    sampling_interval = 1e6/fs
    g = header_geometry(ntraces, dx, shot_number, spi)
    
    spec = segyio.spec()
    spec.ilines  = list(range(ntraces))
    spec.xlines  = [1]
    spec.samples = list(range(nsamples))
    spec.sorting = 2
    spec.format  = 1
    
    with segyio.create(out_sgy, spec) as f:
        
        #populate binary header
        f.bin[3217] = int(sampling_interval)                    #sampling interval in microseconds
        f.bin[3255] = 1                                         #measurement system (1-m / 2-feet)
        
        for i in range(ntraces):
            
            #populate trace headers
            f.header[i][9]      = g['ffid'][i]                  #FFID - original field record number
            f.header[i][17]     = g['spn'][i]                   #energy source point number
            f.header[i][21]     = g['cdp'][i]                   #cdp ensemble number
            f.header[i][71]     = scaler                        #scaler to all coordinates
            f.header[i][73]     = int(100*g['x_src'][i])        #source x coordinate
            f.header[i][77]     = int(100*g['y_src'][i])        #source y coordinate
            f.header[i][81]     = int(100*g['x_rec'][i])        #group x coordinate (receiver)
            f.header[i][85]     = int(100*g['y_rec'][i])        #group y coordinate (receiver)
            f.header[i][109]    = recording_delay               #recording delay
            f.header[i][115]    = nsamples                      #number of samples in this trace
            f.header[i][117]    = int(sampling_interval)        #sample interval in ms for this trace
            f.header[i][181]    = int(100*g['x_cdp'][i])        #CDP X
            f.header[i][185]    = int(100*g['y_cdp'][i])        #CDP Y
            
            #populate the traces with data
            f.trace[i] = np.ascontiguousarray(data[:,i])
            


#%%
def export_sg2(data, fs, dx, shot_number, spi, out_sg2):
    
    print("\nWarning: this function will be depricated, use export_sgy instead")
    
        
    import numpy as np
    import struct
    nsamples, ntraces = np.shape(data)
    
    sampling_interval = 1/fs
    
    #calculate geometry for trace headers
    geometry = header_geometry(ntraces, dx, shot_number, spi)
    
    with open(out_sg2,'wb+') as fid:
    
        # Prepare Strings for File Descriptor Block
        fds1='TRACE_SORT AS_ACQUIRED'                    #File Descriptor String 1
        fds2='UNITS METERS'                              #File Descriptor String 2
    
    
        #Next statement calculates the combined length of strings including the two
        #empty bytes required at the end of the string sub-block to mark the end
        fdslen=len(fds1)+3+len(fds2)+3+2
    
        #Next line checks how many bytes are needed to pad block to a division of 4
        fdspad=4-((fdslen)/4-np.floor((fdslen)/4))/0.25
    
        #Next statement is needed because if tspad==4 then there is no need to pad
        if fdspad==4: 
            fdspad=0 
    
        # Prepare Strings for Trace Descriptor Block
        # This section estimates the lengths of each string in the Trace
        # Descriptor Block and assigns an initial value
    
        tds = {'DELAY ':[],'RECEIVER_LOCATION ':[],'SAMPLE_INTERVAL ':[],'CDP_NUMBER ':[],
            'CDP_TRACE ':[],'RECEIVER_GEOMETRY ':[],'SHOT_SEQUENCE_NUMBER ':[],'SOURCE_LOCATION ':[],
            'SOURCE_GEOMETRY ':[],'SOURCE_STATION_NUMBER ':[],'CHANNEL_NUMBER ':[],'RAW_RECORD ':[],
            'SOURCE VIBROSEIS':['']}
    
        tds['DELAY ']=[("%8.4f"%(geometry['zeros'][0]))]
        tds['RECEIVER_LOCATION ']=[("%15.3f"%(geometry['x_rec'][0]))]
        tds['SAMPLE_INTERVAL ']=[("%0.6f"%(sampling_interval))]
        tds['CDP_NUMBER ']=[("%08d"%(geometry['cdp'][0]))]
        tds['CDP_TRACE ']=[("%08d"%(geometry['x_cdp'][0]))]
        tds['RECEIVER_GEOMETRY ']=[("%15.3f"%(geometry['zeros'][0]))]
        tds['SHOT_SEQUENCE_NUMBER ']=[("%08d"%(geometry['spn'][0]))]
        tds['SOURCE_LOCATION ']=[("%15.3f"%(geometry['x_src'][0]))]
        tds['SOURCE_GEOMETRY ']=[("%15.3f"%(geometry['x_src'][0]))]
        tds['SOURCE_STATION_NUMBER ']=[("%08d"%(geometry['source_number'][0]))]
        tds['CHANNEL_NUMBER ']=[("%08d"%(geometry['channels'][0]))]
        tds['RAW_RECORD ']=[("%08d"%(geometry['ffid'][0]))]
        tds['SOURCE VIBROSEIS']
    
        #Next statement calculates the combined length of strings including the two
        #empty bytes required at the end of the string sub-block to mark the end
    
        lr = [len(list(tds.keys())[ik] + tds[list(tds.keys())[ik]][0]) for ik in range(len(tds))]
    
        tdslen=sum(lr)+3*len(tds)+2
    
        #%Next line checks how many bytes are needed to pad block to a division of 4
        tdspad=4-((tdslen/4)-np.floor((tdslen)/4))/0.25 
    
        #% Next statement is needed because if tspad==4 then there is no need to pad
        if tdspad==4: 
            tdspad=0 
    
    
        #% Write FILE DESCRIPTOR BLOCK
    
        uint16 = '<H'
        uint8 = '<B'
        ulong = '<L'
        float32 = '<f'
    
        fid.write(struct.pack(uint16, int('3A55',base=16)))
        fid.write(struct.pack(uint16, 1))
        fid.write(struct.pack(uint16, ntraces*4))
        fid.write(struct.pack(uint16, ntraces))
        fid.write(struct.pack(uint8, int('01',base=16)))
        fid.write(struct.pack(uint8, 0))
        fid.seek(fid.tell()+1,0)
        fid.write(struct.pack(uint8, int('02',base=16)))
        fid.write(struct.pack(uint8, int('0A',base=16)))
        fid.write(struct.pack(uint8, int('0D',base=16)))
    
    
        #print(fid.read())
        for i in range(1,ntraces+1):
            if i==1: #% When c==1 the pointer skips over the RESERVED Bytes
                fid.seek(fid.tell()+18,0)
                fid.write(struct.pack(ulong, int(32+(4*(ntraces))+fdslen+fdspad+((32+tdslen+tdspad+(nsamples*4))*(i-1)))))
            else:  #%Writes byte location of remaining trace pointer sub-blocks
                fid.write(struct.pack(ulong, int(32+(4*(ntraces))+fdslen+fdspad+((32+tdslen+tdspad+(nsamples*4))*(i-1)))))
    
        # Write Strings for File Descriptor block
        # Write TRACE_SORT string
        fid.write(struct.pack(uint16, int(len(fds1)+3)))
        fid.write(struct.pack(fr'<{len(fds1)}s', bytes(fds1,"utf-8")))
        fid.write(struct.pack(uint8, 0))
        fid.write(struct.pack(uint16, int(len(fds2)+3)))
        fid.write(struct.pack(fr'<{len(fds2)}s', bytes(fds2,"utf-8")))
        fid.write(struct.pack(uint8, 0))
    
        # Write TRACE DESCRIPTOR and DATA BLOCK loops
    
        for i in range(ntraces):
    
            tds['DELAY ']=[("%8.4f"%(geometry['zeros'][i]))]
            tds['RECEIVER_LOCATION ']=[("%15.3f"%(geometry['x_rec'][i]))]
            tds['SAMPLE_INTERVAL ']=[("%0.6f"%(sampling_interval))]
            tds['CDP_NUMBER ']=[("%08d"%(geometry['cdp'][i]))]
            tds['CDP_TRACE ']=[("%08d"%(geometry['x_cdp'][i]))]
            tds['RECEIVER_GEOMETRY ']=[("%15.3f"%(geometry['zeros'][i]))]
            tds['SHOT_SEQUENCE_NUMBER ']=[("%08d"%(geometry['spn'][i]))]
            tds['SOURCE_LOCATION ']=[("%15.3f"%(geometry['x_src'][i]))]
            tds['SOURCE_GEOMETRY ']=[("%15.3f"%(geometry['x_src'][i]))]
            tds['SOURCE_STATION_NUMBER ']=[("%08d"%(geometry['source_number'][i]))]
            tds['CHANNEL_NUMBER ']=[("%08d"%(geometry['channels'][i]))]
            tds['RAW_RECORD ']=[("%08d"%(geometry['ffid'][i]))]
            tds['SOURCE VIBROSEIS']
    
    
            if i==0:
                skip=2+fdspad #Number of bytes to skip to ensure file descriptor
                            #block is a division of 4
            else:
                skip=0 # No skipping is required one first byte of Trace descriptor
                        #block has been written
    
            fid.seek(int(fid.tell()+skip),0)
            fid.write(struct.pack(uint16, int('4422',base=16)))
            fid.write(struct.pack(uint16, int(32+tdslen+tdspad)))
            fid.write(struct.pack(ulong,int(4*nsamples)))
            fid.write(struct.pack(ulong,int(nsamples)))
            fid.write(struct.pack(uint8, int('04',base=16)))
    
    
            for k in range(len(tds)):
    
                par_name = list(tds.keys())[k] + tds[list(tds.keys())[k]][0]
    
                if k==0:
                    skip=19
                else:
                    skip=0
                
                fid.seek(fid.tell()+skip,0)
                fid.write(struct.pack(uint16, int(len(par_name)+3)))
                fid.write(struct.pack(fr'<{len(par_name)}s', bytes(par_name,"utf-8")))
                fid.write(struct.pack(uint8, 0))
    
    
            for j in range(nsamples):
                if j==0:
                    fid.seek(int(fid.tell()+2+tdspad),0)
                    fid.write(struct.pack(float32, data[j,i]))
                else:
                    fid.write(struct.pack(float32, data[j,i]))

    fid.close()
