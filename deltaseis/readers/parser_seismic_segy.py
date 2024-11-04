# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:49:08 2023
Tip: in powershell use: format-hex <segy_file_path> > binairy.txt to view
the segy file in binairy format. Scroll to the end to look how many junk bytes 
there are.

@author: nieboer
"""

def parse_segy(segy_filepath):
    ''' Deduces how the segy file should be opened: 
        (1) first trying with big endian byte order, if this fails
        (2) secondly tries with little endian byte order, if this fails
        (3) tries to fix the input segy by removing null bytes at the end
        of a copy of the file with _TEMPORARY appended to the file name. If
        this last option is required, the class will be instantiated with
        the _TEMPORARY file version.
        
        Parameters
        ----------
        segy_filepath: String or pathlib.Path object
        
        Returns
        -------
        segy_filepath_out: String or pathlib.Path object
        Returns the filepath to be used which is the original filepath or if
        null bytes were remove, the filepath _CLEANED appended to it
        
        endian: String
        the byte order of the segy file 'big' or 'little'
        '''
        
    from pathlib import Path
    import struct
    segy_filepath = Path(segy_filepath)
        
    try: 
        endian, byte_size, data_sample_format = get_endianness(segy_filepath)
        print('The segy file has {} endian byte order'.format(endian))
        segy_filepath_out = segy_filepath
          
    except NameError:
        segy_filepath_out = segy_filepath.with_stem(segy_filepath.stem + "_TEMPORARY")
                
        with open(segy_filepath, 'rb') as file_in:
            with open(segy_filepath_out, 'wb') as file_out:
                data = file_in.read()
                file_byte_length = len(data)
                
                #trick to check endianness from bin header
                file_in.seek(3224, 0)
                format_code_bin = file_in.read(2)
                format_code = struct.unpack('>h', format_code_bin)[0]
                
                if 1 <= format_code <= 4:
                    endian = 'big'
                    mode = '>h'
                else:
                    endian = 'little'
                    mode = '<h'
        
                    
                #read number of samples and sample format from the binary header               
                file_in.seek(3220, 0)
                number_of_samples_bin = file_in.read(2)
                number_of_samples = struct.unpack(mode, number_of_samples_bin)[0]
                             
                file_in.seek(3224, 0)
                format_code_bin = file_in.read(2)
                format_code = struct.unpack(mode, format_code_bin)[0]
                
                #retrieve sample byte size from format_code 
                byte_size, data_sample_format = get_data_sample_format(format_code)
                
               #determine the number of traces that fit in the total file size
                text_header_bytes   = 3200
                bin_header_bytes    = 400
                trace_header_bytes  = 240
                
                full_trace_byte_length = trace_header_bytes + byte_size*number_of_samples
                number_of_traces = int((len(data) - text_header_bytes - bin_header_bytes) / full_trace_byte_length)
                
                #now determin the amount of junk bytes at the end of the file
                number_of_junk_bytes = len(data) - text_header_bytes - bin_header_bytes - number_of_traces*full_trace_byte_length
                
                #remove the junk bytes from the end of the file        
                data = data[:-number_of_junk_bytes]
                    
                file_out.write(data)                   
                print("{} null bytes removed from end of file".format(file_byte_length - len(data)))
                print('The cleaned segy file has {} endian byte order and has sample format {}'.format(endian, data_sample_format))
                
    return segy_filepath_out, endian, data_sample_format

def get_data_sample_format(format_code):
    'retrieves the sample format based on the seg-y format code'
    
    #retrieve sample byte size from format_code ### make this separate def
    if format_code==1 or format_code==5:
        byte_size = 4
        data_sample_format = 'float32'
    elif format_code==2:
        byte_size = 4
        data_sample_format = 'int32'
    elif format_code==3:
        byte_size = 2
        data_sample_format = 'int16'
    elif format_code==6:
        byte_size = 8
        data_sample_format = 'float64'
    elif format_code==9:
        byte_size = 8
        data_sample_format = 'int64'
        
    return byte_size, data_sample_format
                    
def get_endianness(segy_filepath):
    '''determine if the input segy file has big or little endian byte order
    or print error if it is neither'''
    
    import segyio
    
    try: 
        with segyio.open(segy_filepath, ignore_geometry=True, endian='big') as big:
            endian= 'big'
            byte_size, data_sample_format = get_data_sample_format(big.bin[3225])
        
    except RuntimeError:
        try: 
            with segyio.open(segy_filepath, ignore_geometry=True, endian='little') as little:
                endian= 'little'
                byte_size, data_sample_format = get_data_sample_format(little.bin[3225])
        except RuntimeError:
            print("Segy file could not be recognized as either having 'big' or 'little' endian byte order")
            
    return endian, byte_size, data_sample_format