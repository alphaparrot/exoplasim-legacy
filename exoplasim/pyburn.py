"""
Read raw exoplasim output files and postprocess them into netCDF output files.
"""
import numpy as np
import netCDF4 as nc
import struct

'''
This module is intended to be a near-replacement for the C++ burn7 utility, which in its present
form relies on deprecated and unsupported versions of the netcdf C and C++ libraries. This utility
relies instead on the netCDF4-python package, whose usage pattern is unlikely to become deprecated
at any point in the near future. Some features of burn7 are not included; notably at present only
netCDF output is supported, rather than the additional SRV and GRaDs formats supported by burn7.
Additionally, at least at first, there are the following restrictions:

    **Vertical levels will only be output in native sigma levels,** and will not be interpolated onto
    a pressure grid. This is unlikely to be supported in the future either, as the conversion of sigma
    to pressure is trivial, and interpolating onto a different pressure grid at the postprocessing step
    risks misleading end users and is liable to run into problems with non-Earthlike surface pressures.
    
    **Horizontal grid will only be Gaussian-spaced latitude-longitude,** rather than being optionally
    spherical harmonics, fourier coefficients, or zonal averages. The additional horizontal options present
    in burn7 will be supported in pyburn in future versions.
    
    **The Mars shortcut for derived variables is unsupported.** This will be supported in future versions.
'''

def _getEndian(fbuffer):
    '''Determine Endian-ness of the buffer'''
    import sys
    byteorder = sys.byteorder
    if byteorder=="little":
        endian="<"
    else:
        endian=">"
    tag = struct.unpack('i',fbuffer[:4])[0]
    if tag>1024:   #If we get a huge number for what should be an 8-word header, flip the byte order
        if endian==">": #We're in a Big-Endian system, and the buffer was written in little-Endian.
            endian="<"
        else: #We're in a little-Endian system, and the buffer was written in Big-Endian.
            endian=">"
    return endian

def _getmarkerlength(fbuffer,en):
    '''Determine how many bytes were used for record markers'''
    tag = struct.unpack(en+'i',fbuffer[:4])[0]
    markerlength=4
    if tag==0:
        markerlength=8
    return markerlength

def _getwordlength(fbuffer,n,en,fmt='i'):
    '''Determine if we're dealing with 32-bit output or 64-bit.
    
    At the moment, all exoplasim outputs are 32-bit, even if run in 64-bit. However, we should build
    compatibility for this eventuality in the future.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with `mode='rb'` and read with `file.read()`. 
    n : int
        The index of the byte at which to check. This should be the start of the first word of the 
        variable in question.
    en : str
        Endianness, denoted by ">" or "<"
    fmt : str
        Expected type of the word--'i' for integers (longs), and 'd' for floats (doubles). 
    
    Returns
    -------
    int, str
       Word length in bytes, and format string for a word--4 for 32 bit, and 8 for 64 bit. 'i' for a 
       4-byte int, 'l' for an 8-byte int, 'f' for a 4-byte float, and 'd' for an 8-byte float.
    '''

    tag = struct.unpack(en+fmt,fbuffer[n:n+4])[0]
    wordlength=4
    if tag==0:
        wordlength=8
    if fmt=='i' and wordlength==8:
        fmt='l'
    elif fmt=='f' and wordlength==8:
        fmt='d'
    return wordlength,fmt

def _getknownwordlength(fbuffer,n,en,ml):
    '''Determine word length of a data variable in cases where we know that the header is 8 words and is 32-bit.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with `mode='rb'` and read with `file.read()`. 
    n : int
        The index of the word at which to start, in bytes. A 32-bit word has length 4, so the current 
        position in words would be 4*n assuming 4-byte words, or 8*n if 64 bits and 8-byte words.
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
    
    Returns
    -------
    int, str
       Word length in bytes, and format string for a word--4 for 32 bit, and 8 for 64 bit. 
       'f' for a 4-byte float, and 'd' for an 8-byte float.
    '''
    
    htag = struct.unpack(en+'i',fbuffer[n:n+ml])
    n+=ml
    header = struct.unpack(en+8*'i',fbuffer[n:n+32])
    n+=32+ml #Add one word for restatement of header length
    dtag = struct.unpack(en+'i',fbuffer[n:n+ml])[0]
    
    dim1 = header[4]
    dim2 = header[5]
    
    length = dim1*dim2
    
    wordlength = dtag//length #dtag tells us the length of the coming record in bytes, while length is the
                              #length of the coming record in words. The ratio of the two is thus wordlength.
    if wordlength==4:
        fmt='f'
    else:
        fmt='d'
    return wordlength,fmt

def readrecord(fbuffer,n,en,ml):
    '''Read a Fortran record from the buffer, starting at index n, and return the header, data, and updated n.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with `mode='rb'` and read with `file.read()`. 
    n : int
        The index of the word at which to start, in bytes. A 32-bit word has length 4, so the current 
        position in words would be 4*n assuming 4-byte words, or 8*n if 64 bits and 8-byte words.
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
        
    Returns
    -------
    array-like, array-like, int
        A tuple containing first the header, then the data contained in the record, and finally the new
        position in the buffer in bytes.
    '''
    if n<len(fbuffer):
        wl,fmt = _getknownwordlength(fbuffer,n,en,ml)
                
        headerlength = int(struct.unpack(en+'i',fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4)])
        n+=headerlength*4+ml #Add one word for restatement of header length (for backwards seeking)
        datalength = int(struct.unpack(en+'i',fbuffer[n:n+ml])[0]//wl)
        n+=ml
        data = struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl])
        n+=datalength*wl+ml #additional 4 for restatement of datalength
        return header,data,n
    else:
        raise Exception("Reached end of buffer!!!")
    
def readvariablecode(fbuffer,kcode,en,ml):
    '''Seek through a binary output buffer and extract all records associated with a variable code.
    
    Note, assembling a variable list piece by piece in this way may be slower than reading **all** variables
    at once, because it requires seeking all the way through the buffer multiple times for each variable.
    This will likely only be faster if you only need a small number of variables.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with `mode='rb'` and read with `file.read()`.
    kcode : int
        The integer code associated with the variable. For possible codes, refer to the 
        `Postprocessor Variable Codes. <postprocessor.html#postprocessor-variable-codes>`_
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
    
    Returns
    -------
    array-like, array-like
        A tuple containing first the header, then the variable data, as one concatenated 1D variable.
    '''
    n = 0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml)
    
    variable = None
    
    while n<len(fbuffer):
        
        recordn0 = n
        
        headerlength = int(struct.unpack(en+'i',fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4)])
        n+=headerlength*4+ml
        datalength = struct.unpack(en+'i',fbuffer[n:n+ml])[0]
        n+=ml
        if header[0]==kcode:
            dataheader = header
            wl, fmt = _getknownwordlength(fbuffer,recordn0,en,ml)
            datalength = int(datalength//wl)
            if not variable:
                variable = np.array(struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            else:
                variable = np.append(variable,struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            n+=datalength*wl+ml
        else: #Fast-forward past this variable without reading it.
            n+=datalength+ml
    
    return dataheader, variable
    
def readallvariables(fbuffer):
    '''Extract all variables and their headers from a file byte buffer.
    
    Doing this and then only keeping the codes you want may be faster than extracting variables one by one,
    because it only needs to seek through the file one time.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with `mode='rb'` and read with `file.read()`.
    
    Returns
    -------
    dict, dict
        A dictionary containing all variable headers (by variable code), and a dictionary containing all
        variables, again by variable code.
    '''
    
    en = _getEndian(fbuffer)
    ml = _getmarkerlength(fbuffer,en)
    
    n=0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml)
    
    headers= {'main':mainheader}
    variables = {'main':zsig}
    nlev=mainheader[6]
    variables["sigmah"] = zsig[:nlev]
    
    while n<len(fbuffer):
        header,field,n = readrecord(fbuffer,n,en,ml)
        kcode = str(header[0])
        if kcode not in variables:
            variables[kcode] = np.array(field)
            headers[kcode] = header
            print("variable %d"%kcode)
        else:
            variables[kcode] = np.append(variables[kcode],field)
    
    return headers, variables
    
    
def refactorvariable(variable,header,nlev=10):
    '''Given a 1D data array extracted from a file with :py:func:`readrecord <exoplasim.pyburn.readrecord>`, reshape it into its appropriate dimensions.
    
    Parameters
    ----------
    variable : array-like
        Data array extracted from an output file using :py:func:`readrecord <exoplasim.pyburn.readrecord>`.
        Can also be the product of a concatenated file assembled with
        :py:func:`readvariable <exoplasim.pyburn.readvariable>`.
    header : array-like
        The header array extracted from the record associated with `variable`. This header contains
        dimensional information.
    nlev : int, optional
        The number of vertical levels in the variable. If 1, vertical levels will not be a dimension in
        the output variable.
        
    Returns
    -------
    numpy.ndarray
        A numpy array with dimensions (time,lat,lon) if `nlevs=1`, or (time,lev,lat,lon) otherwise.
    '''
    if header[1]==1:
        nlevs=nlev
    else:
        nlevs=1
    dim1 = header[4]
    dim2 = header[5]
    ntimes = int(len(variable)//(dim1*dim2*nlevs))
    if nlevs==1:
        newvar = np.reshape(variable,(ntimes,dim2,dim1))
    else:
        newvar = np.reshape(variable,(ntimes,nlevs,dim2,dim1))
    return newvar

def readfile(filename):
    '''Extract all variables from a raw plasim output file and refactor them into the right shapes
    
    Parameters
    ----------
    filename : str
        Path to the output file to read
        
    Returns
    -------
    dict
        Dictionary of model variables, indexed by numerical code
    '''
    
    with open(filename,"rb") as fb:
        fbuffer = fb.read()
    
    headers, variables = readallvariables(fbuffer)
    
    nlevs = len(variables['sigmah'])
    
    kcodes = list(variables.keys())
    kcodes.remove('sigmah')
    
    data = {}
    
    for key in kcodes:
        data[key] = refactorvariable(variables[key],headers[key],nlev=nlevs)
        
    # Add in latitude-longitude generation, plus mid-level computation, and time array    
        
    return data
