"""
Read raw exoplasim output files and postprocess them into netCDF output files.
"""
import numpy as np
import netCDF4 as nc
import struct
import pyfft
import gcmt
import scipy

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

#dictionary that can be searched by string (integer) codes
ilibrary = {"110":["mld"  ,"mixed_layer_depth"               ,"m"          ], 
            "129":["sg"   ,"surface_geopotential"            ,"m2 s-2"     ],
            "130":["ta"   ,"air_temperature"                 ,"K"          ],
            "131":["ua"   ,"eastward_wind"                   ,"m s-1"      ],
            "132":["va"   ,"northward_wind"                  ,"m s-1"      ],
            "133":["hus"  ,"specific_humidity"               ,"1"          ],
            "134":["ps"   ,"surface_air_pressure"            ,"hPa"        ],
            "135":["wap"  ,"vertical_air_velocity"           ,"Pa s-1"     ],
            "137":["wa"   ,"upward_wind"                     ,"m s-1"      ],
            "138":["zeta" ,"atm_relative_vorticity"          ,"s-1"        ],
            "139":["ts"   ,"surface_temperature"             ,"K"          ],
            "140":["mrso" ,"lwe_of_soil_moisture_content"    ,"m"          ],
            "141":["snd"  ,"surface_snow_thickness"          ,"m"          ],
            "142":["prl"  ,"lwe_of_large_scale_precipitation","m s-1"      ],
            "143":["prc"  ,"convective_precipitation_rate"   ,"m s-1"      ],
            "144":["prsn" ,"lwe_of_snowfall_amount"          ,"m s-1"      ],
            "145":["bld"  ,"dissipation_in_atmosphere_bl"    ,"W m-2"      ],
            "146":["hfss" ,"surface_sensible_heat_flux"      ,"W m-2"      ],
            "147":["hfls" ,"surface_latent_heat_flux"        ,"W m-2"      ],
            "148":["stf"  ,"streamfunction"                  ,"m2 s-2"     ],
            "149":["psi"  ,"velocity_potential"              ,"m2 s-2"     ],
            "151":["psl"  ,"air_pressure_at_sea_level"       ,"hPa"        ],
            "152":["pl"   ,"log_surface_pressure"            ,"1"          ],
            "155":["d"    ,"divergence_of_wind"              ,"s-1"        ],
            "156":["zg"   ,"geopotential_height"             ,"m"          ],
            "157":["hur"  ,"relative_humidity"               ,"1"          ],
            "158":["tps"  ,"tendency_of_surface_air_pressure","Pa s-1"     ],
            "159":["u3"   ,"ustar"                           ,"m3 s-3"     ],
            "160":["mrro" ,"surface_runoff"                  ,"m s-1"      ],
            "161":["clw"  ,"liquid_water_content"            ,"1"          ],
            "162":["cl"   ,"cloud_area_fraction_in_layer"    ,"1"          ],
            "163":["tcc"  ,"total_cloud_cover"               ,"1"          ],
            "164":["clt"  ,"cloud_area_fraction"             ,"1"          ],
            "165":["uas"  ,"eastward_wind_10m"               ,"m s-1"      ],
            "166":["vas"  ,"northward_wind_10m"              ,"m s-1"      ],
            "167":["tas"  ,"air_temperature_2m"              ,"K"          ],
            "168":["td2m" ,"dew_point_temperature_2m"        ,"K"          ],
            "169":["tsa"  ,"surface_temperature_accumulated" ,"K"          ],
            "170":["tsod" ,"deep_soil_temperature"           ,"K"          ],
            "171":["dsw"  ,"deep_soil_wetness"               ,"1"          ],
            "172":["lsm"  ,"land_binary_mask"                ,"1"          ],
            "173":["z0"   ,"surface_roughness_length"        ,"m"          ],
            "174":["alb"  ,"surface_albedo"                  ,"1"          ],
            "175":["as"   ,"surface_albedo"                  ,"1"          ],
            "176":["rss"  ,"surface_net_shortwave_flux"      ,"W m-2"      ],
            "177":["rls"  ,"surface_net_longwave_flux"       ,"W m-2"      ],
            "178":["rst"  ,"toa_net_shortwave_flux"          ,"W m-2"      ],
            "179":["rlut" ,"toa_net_longwave_flux"           ,"W m-2"      ],
            "180":["tauu" ,"surface_eastward_stress"         ,"Pa"         ],
            "181":["tauv" ,"surface_northward_stress"        ,"Pa"         ],
            "182":["evap" ,"lwe_of_water_evaporation"        ,"m s-1"      ],
            "183":["tso"  ,"climate_deep_soil_temperature"   ,"K"          ],
            "184":["wsoi" ,"climate_deep_soil_wetness"       ,"1"          ],
            "199":["vegc" ,"vegetation_cover"                ,"1"          ],
            "203":["rsut" ,"toa_outgoing_shortwave_flux"     ,"W m-2"      ],
            "204":["ssru" ,"surface_solar_radiation_upward"  ,"W m-2"      ],
            "205":["stru" ,"surface_thermal_radiation_upward","W m-2"      ],
            "207":["tso2" ,"soil_temperature_level_2"        ,"K"          ],
            "208":["tso3" ,"soil_temperature_level_3"        ,"K"          ],
            "209":["tso4" ,"soil_temperature_level_4"        ,"K"          ],
            "210":["sic"  ,"sea_ice_cover"                   ,"1"          ],
            "211":["sit"  ,"sea_ice_thickness"               ,"m"          ],
            "212":["vegf" ,"forest_cover"                    ,"1"          ],
            "218":["snm"  ,"snow_melt"                       ,"m s-1"      ],
            "221":["sndc" ,"snow_depth_change"               ,"m s-1"      ],
            "230":["prw"  ,"atmosphere_water_vapor_content"  ,"kg m-2"     ],
            "232":["glac" ,"glacier_cover"                   ,"1"          ],
            "238":["tsn"  ,"snow_temperature"                ,"K"          ],
            "259":["spd"  ,"wind_speed"                      ,"m s-1"      ],
            "260":["pr"   ,"total_precipitation"             ,"m s-1"      ],
            "261":["ntr"  ,"net_top_radiation"               ,"W m-2"      ],
            "262":["nbr"  ,"net_bottom_radiation"            ,"W m-2"      ],
            "263":["hfns" ,"surface_downward_heat_flux"      ,"W m-2"      ],
            "264":["wfn"  ,"net_water_flux"                  ,"m s-1"      ],
            "265":["dqo3" ,"ozone_concentration"             ,"unknown"    ],
            "266":["lwth" ,"local_weathering"                ,"W_earth"    ],
            "267":["grnz" ,"ground_geopotential"             ,"m2 s-2"     ],
            "273":["dpdx" ,"d(ps)/dx"                        ,"Pa m-1"     ],
            "274":["dpdy" ,"d(ps)/dy"                        ,"Pa m-1"     ],
            "277":["hlpr" ,"half_level_pressure"             ,"Pa"         ],
            "278":["flpr" ,"full_level_pressure"             ,"Pa"         ],
            "301":["icez" ,"glacier_geopotential"            ,"m2 s-2"     ],
            "302":["netz" ,"net_geopotential"                ,"m2 s-2"     ],
            "318":["czen" ,"cosine_solar_zenith_angle"       ,"nondimen"   ],
            "319":["wthpr","weatherable_precipitation"       ,"mm day-1"   ],
            "320":["mint" ,"minimum_temperature"             ,"K"          ],
            "321":["maxt" ,"maximum_temperature"             ,"K"          ],
            "322":["cape" ,"convective_available_potential_energy","J kg-1"],
            "323":["lnb"  ,"level_neutral_buoyancy"          ,"hPa"        ],
            "324":["sdef" ,"troposphere_entropy_deficit"     ,"1"          ],
            "325":["absz" ,"sigma=0.85_abs_vorticity"        ,"s-1"        ],
            "326":["umax" ,"maximum_potential_intensity"     ,"m s-1"      ],
            "327":["vent" ,"ventilation_index"               ,"1"          ],
            "328":["vrumax","ventilation_reduced_maximum_wind","m s-1"     ],
            "329":["gpi"  ,"genesis_potential_index"            ,"1"       ],
            "404":["dfu"  ,"shortwave_up"                    ,"W m-2"      ],
            "405":["dfd"  ,"shortwave_down"                  ,"W m-2"      ],
            "406":["dftu"  ,"longwave_up"                    ,"W m-2"      ],
            "407":["dftd"  ,"longwave_down"                  ,"W m-2"      ],
            "408":["dtdt"  ,"rad_heating_rate"               ,"K s-1"      ],
            "409":["dfdz"  ,"flux_convergence"               ,"W m-3"      ]}

#dictionary that can be searched by string codes
slibrary = {}
for key in ilibrary:
    kcode = int(key)
    name = ilibrary[key][0]
    longn = ilibrary[key][1]
    units = ilibrary[key][2]
    slibrary[name] = [kcode,longn,units]
                    
geopotcode  = 129
tempcode    = 130
ucode       = 131
vcode       = 132
humcode     = 133
pscode      = 134
omegacode   = 135
wcode       = 137
vortcode    = 138
tscode      = 139
stfcode     = 148
vpotcode    = 149
slpcode     = 151
lnpscode    = 152
divcode     = 155
geopotzcode = 156
rhumcode    = 157
spdcode     = 259
preccode    = 260
ntopcode    = 261
nbotcode    = 262
nheatcode   = 263
nh2ocode    = 264
swatmcode   = 268
lwatmcode   = 269
natmcode    = 270
sruncode    = 271
dpsdxcode   = 273
dpsdycode   = 274
freshcode   = 275
            
hpresscode  = 277
fpresscode  = 278
thetahcode  = 279
thetafcode  = 280

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
    elif tag==0: #it was an 8-byte record marker
        tag = struct.unpack('l',fbuffer[:8])[0]
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

def _getknownwordlength(fbuffer,n,en,ml,mf):
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
    mf : str
        Format of the record marker ('i' or 'l')
    
    Returns
    -------
    int, str
       Word length in bytes, and format string for a word--4 for 32 bit, and 8 for 64 bit. 
       'f' for a 4-byte float, and 'd' for an 8-byte float.
    '''
    
    htag = struct.unpack(en+mf,fbuffer[n:n+ml])
    n+=ml
    header = struct.unpack(en+8*'i',fbuffer[n:n+32])
    n+=32+ml #Add one word for restatement of header length
    dtag = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
    
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

def _ql(k,p):
    '''Part of gaussian latitude computation. Translated from exoplasim fortran.'''
    z0 = np.arccos(p)
    z1 = 1.0
    z2 = 0.0
    for j in range(k , -1 , -2):
        z3 = z1 * np.cos(z0 * j)
        z2 = z2 + z3
        z4 = (k-j+1) * (k+j) * 0.5
        z1 = z1 * z4 / (z4 + (j-1))
    if (k%2 == 0):
        z2 = z2 - 0.5 * z3

    z0 = np.sqrt(2.0)
    for j in range(1 ,k+1):
        z0 = z0 * np.sqrt(1.0 - 0.25 / (j*j))
    return z0 * z2

def _qld(k,p):
    '''Part of gaussian latitude computation. Translated from exoplasim fortran.'''
    z = p * _ql(k,p) - np.sqrt((k + k + 1.0) / (k + k - 1.0)) * _ql(k-1,p)
    return (p * p - 1.0) / (k * z)

def _gaulat(klat):
    '''Given klat latitudes, compute the Gaussian-spaced latitudes and return sin(lat) and the gaussian weights.
    
    Parameters
    ----------
    klat : int
        Number of latitudes in the model
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Sine of the latitudes, and the associated gaussian weights.
    '''
    NITER = 50
    ZEPS = 1.0e-16
    z0 = np.pi / (2*klat+1)
    z1 = 1.0 / (klat*klat*8)
    z4 = 2.0 / (klat*klat)

    pz0 = np.zeros(klat)
    pzw = np.zeros(klat)

    for jlat in range(1 , klat//2+1):
        z2 = z0 * (2*jlat - 0.5)
        z2 = np.cos(z2 + z1 / np.tan(z2))
        for jiter in range(0 , NITER):
            z3 = _ql(klat,z2) * _qld(klat,z2)
            z2 = z2 - z3
            if (abs(z3) < ZEPS):
                break
        z5 = _ql(klat-1,z2) / np.sqrt(klat - 0.5)
        pz0[jlat-1] = z2
        pzw[jlat-1] = z4 * (1.0 - z2 * z2) / (z5 * z5)
        pz0[klat-1-jlat+1] = -z2
        pzw[klat-1-jlat+1] = pzw[jlat-1]
    return pz0,pzw #sid,gwd


def _legini(sid,gwd,NLAT,NTRU):
    '''Compute coefficients for Legendre transformations
    
    Parameters
    ----------
    sid : array-like
        Sine of latitudes
    gwd : array-like
        Gaussian latitude weights
    NTRU : int
        Truncation number
    '''
    NLAT = len(sid)
    NLPP = NLAT
    NCSP = (NTRU+1)*(NTRU+2)//2
    qi = np.zeros((NLPP,NCSP)) # ! P(m,n) = Associated Legendre Polynomials
    qj = np.zeros((NLPP,NCSP)) # ! Q(m,n) = Used for d/d(mu)
    qc = np.zeros((NLPP,NCSP)) # ! P(m,n) * gwd              used in fc2sp
    qu = np.zeros((NLPP,NCSP)) # ! Q(mn,) * gwd / cos2       used in mktend
    qv = np.zeros((NLPP,NCSP)) # ! P(m,n) * gwd / cos2 * m   used in mktend
    qe = np.zeros((NLPP,NCSP)) # ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
    qq = np.zeros((NLPP,NCSP)) # ! P(m,n) / (n*(n+1)) * m    used in dv2uv
    qm = np.zeros((NLPP,NCSP)) # ! Q(m,n) / (n*(n+1))        used in dv2uv
    zpli = np.zeros(NCSP)
    zpld = np.zeros(NCSP)
    for jlat in range(1 , NLPP+1):

#! set p(0,0) and p(0,1)

        zgwd    = gwd[jlat-1]#            ! gaussian weight - from inigau
        zsin    = sid[jlat-1]#            ! sin(phi) - from inigau
        zcsq    = 1.0 - zsin * zsin  #! cos(phi) squared
        zgwdcsq = zgwd / zcsq  #        ! weight / cos squared
        f1m     = np.sqrt(1.5)
        zpli[0] = np.sqrt(0.5)
        zpli[1] = f1m * zsin
        zpld[0] = 0.0
        lm      = 2

        #! loop over wavenumbers

        for m in range(0 , NTRU+1):
            if (m > 0):
                lm  = lm + 1
                f2m = -f1m * np.sqrt(zcsq / (m+m))
                f1m =  f2m * np.sqrt(m+m + 3.0)
                zpli[lm-1] = f2m
                if (lm < NCSP):
                    lm = lm + 1
                    zpli[lm-1] =       f1m * zsin
                    zpld[lm-2] =  -m * f2m * zsin
                #endif ! (lm < NCSP)
            #endif ! (m > 0)

            amsq = m * m

            for n in range(m+2 , NTRU+1):
                lm = lm + 1
                z1 = np.sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
                z2 = zsin * zpli[lm-2] - z1 * zpli[lm-3]
                zpli[lm-1] = z2 * np.sqrt((4*n*n-1) / (n*n-amsq))
                zpld[lm-2] = (1-n) * z2 + n * z1 * zpli[lm-3]
            #enddo ! n

            if (lm < NCSP): #then ! mode (m,NTRU)
                z3 = np.sqrt((NTRU*NTRU-amsq) / (4*NTRU*NTRU-1))
                zpld[lm-1]=-NTRU*zsin*zpli[lm-1] + (NTRU+NTRU+1)*zpli[lm-2]*z3
            else: #               ! mode (NTRU,NTRU)
                zpld[lm-1]=-NTRU*zsin*zpli[lm-1]
            #endif
        #enddo ! m

        lm = 0
        for m in range(0 , NTRU+1):
            for n in range(m , NTRU+1):
                lm = lm + 1
                znn1 = 0.0
                if (n > 0):
                    znn1 = 1.0 / (n*(n+1))
                qi[jlat-1,lm-1] = zpli[lm-1]
                qj[jlat-1,lm-1] = zpld[lm-1]
                qc[jlat-1,lm-1] = zpli[lm-1] * zgwd
                qu[jlat-1,lm-1] = zpli[lm-1] * znn1 * m
                qv[jlat-1,lm-1] = zpld[lm-1] * znn1
                qe[jlat-1,lm-1] = zpld[lm-1] * zgwdcsq
                qq[jlat-1,lm-1] = zpli[lm-1] * zgwdcsq * n * (n+1) * 0.5
                qm[jlat-1,lm-1] = zpli[lm-1] * zgwdcsq * m
    
    return qi,qj,qc,qu,qv,qe,qq,qm

    
def _sp2fc(sp,NLON,NLAT,NTRU,qi):
    '''Convert spectral harmonics into fourier coefficients
    
    Parameters
    ----------
    sp : array-like
        Spectral field
    NLON : int 
        number of longitudes
    NLAT : int
        number of latitudes
    NTRU : int
        truncation number
    qi : array-like
        Legendre polynomials
    '''
    NTP1 = NTRU + 1
    NLPP = NLAT
    fc = np.zeros((NLPP,NLON//2,2))               #(2,NLON/2,NLPP) ! Fourier coefficients
    NCSP = (NTRU+1)*(NTRU+2)//2
    sp = np.reshape(sp,(NCSP,2)) #Presumably we had NCSP*2 length previously....?
#integer :: l ! Loop index for latitude
#integer :: m ! Loop index for zonal wavenumber m
#integer :: n ! Loop index for total wavenumber n
#integer :: w ! Loop index for spectral mode

    for l in range(1 , NLPP+1):
        w = 1  
        for m in range(1 , NTP1+1):
            for n in range(m , NTP1+1):
                fc[l-1,m-1,0] = fc[l-1,m-1,0] + qi[l-1,w-1] * sp[w-1,0]#!*skspgp(n)
                fc[l-1,m-1,1] = fc[l-1,m-1,1] + qi[l-1,w-1] * sp[w-1,1]#!*skspgp(n)
                w = w + 1
    return fc


def _sp3fc(sp,NLON,NLAT,NTRU,NLEV,qi):
    '''Convert spectral harmonics into fourier coefficients for a 3D field
    
    Parameters
    ----------
    sp : array-like
        Spectral field
    NLON : int 
        number of longitudes
    NLAT : int
        number of latitudes
    NTRU : int
        truncation number
    NLEV : int 
        Number of vertical levels
    qi : array-like
        Legendre polynomials
    '''
    fc = np.zeros((NLEV,NLAT,NLON//2,2))
    for jlev in range(NLEV):
        fc[jlev,...] = _sp2fc(sp[jlev,...],NLON,NLAT,NTRU,qi)
    return fc



def _ifft8(fc,n,lot):
    SQRT2 = 1.414213562373095
    #dimension a(n*lot),c(n*lot)
    a = np.reshape(fc,(n*lot))
    la = n // 8

    c = np.zeros(n*lot)
    
    for i in range(0,la*lot):
        i0 = (i//la) * n + i%la + 1
        i1 = i0 + la
        i2 = i1 + la
        i3 = i2 + la
        i4 = i3 + la
        i5 = i4 + la
        i6 = i5 + la
        i7 = i6 + la

        a0p7 = a[i0-1] + a[i7-1]
        a0m7 = a[i0-1] - a[i7-1]
        a1p5 = a[i1-1] + a[i5-1]
        a1m5 = a[i1-1] - a[i5-1]
        a2p6 = a[i2-1] + a[i6-1]
        a2m6 = a[i2-1] - a[i6-1]

        a0p7p3   = a0p7 + a[i3-1]
        a0p7m3   = a0p7 - a[i3-1]
        a0m7p4   = 2.0 * (a0m7 + a[i4-1])
        a0m7m4   = 2.0 * (a0m7 - a[i4-1])
        a1m5p2p6 = SQRT2 * (a1m5 + a2p6)
        a1m5m2p6 = SQRT2 * (a1m5 - a2p6)

        c[i0-1]  = 2.0 * (a0p7p3 + a1p5)
        c[i2-1]  = 2.0 * (a0p7m3 - a2m6)
        c[i4-1]  = 2.0 * (a0p7p3 - a1p5)
        c[i6-1]  = 2.0 * (a0p7m3 + a2m6)

        c[i1-1]  = a0m7m4 + a1m5m2p6
        c[i3-1]  = a0m7p4 - a1m5p2p6
        c[i5-1]  = a0m7m4 - a1m5m2p6
        c[i7-1]  = a0m7p4 + a1m5p2p6
    return np.reshape(c,(lot,n))

def _ifft4(fc,trigs,n,lot,la):
    
    a = np.zeros((lot,n))
    c = np.copy(fc)
    
    if (la == 1):
        a[:,0] = 0.5 * c[:,0]
        a[:,n-1] = 0.0
        a[:,1:n-1] = c[:,2:n]
    else:
        a = c

    kstop=(n-4)//8

    i1 = n//2 - la
    i2 = n   - la
    i5 = i1  + la

    j1 = la
    j2 = la+j1
    j3 = la+j2
    j4 = la+j3
    j5 = la+j4
    j6 = la+j5
    j7 = la+j6

    for i in range(1,la+1):
        for l in range(1,lot+1):
            c[l-1,   i-1] = a[l-1,i-1] + a[l-1,i2+i-1] + a[l-1,i1+i-1]
            c[l-1,j1+i-1] = a[l-1,i-1] - a[l-1,i2+i-1] - a[l-1,i5+i-1]
            c[l-1,j2+i-1] = a[l-1,i-1] + a[l-1,i2+i-1] - a[l-1,i1+i-1]
            c[l-1,j3+i-1] = a[l-1,i-1] - a[l-1,i2+i-1] + a[l-1,i5+i-1]

    iink  = 2 * la
    jbase = 4 * la + 1
    i0    = la
    i1    = i0 + n//2
    i2    = n - 3 * la
    i3    = i2 - n//2
    i4    = i0 + la
    i5    = i1 + la
    i6    = i2 + la
    i7    = i3 + la

    for k in range(la,kstop+1,la):
        kb=k+k
        kc=kb+kb
        kd=kc+kb
        c1=trigs[kb  ]
        s1=trigs[kb+1]
        c2=trigs[kc  ]
        s2=trigs[kc+1]
        c3=trigs[kd  ]
        s3=trigs[kd+1]
        for i in range(1 , la+1):
            j = jbase
            for l in range(1,lot+1):
                a0p2 = a[l-1,i0+i-1] + a[l-1,i2+i-1]
                a0m2 = a[l-1,i0+i-1] - a[l-1,i2+i-1]
                a1p3 = a[l-1,i1+i-1] + a[l-1,i3+i-1]
                a1m3 = a[l-1,i1+i-1] - a[l-1,i3+i-1]
                a4p6 = a[l-1,i4+i-1] + a[l-1,i6+i-1]
                a4m6 = a[l-1,i4+i-1] - a[l-1,i6+i-1]
                a5p7 = a[l-1,i5+i-1] + a[l-1,i7+i-1]
                a5m7 = a[l-1,i5+i-1] - a[l-1,i7+i-1]

                a0p2m1p3 = a0p2 - a1p3
                a4m6m5m7 = a4m6 - a5m7

                c[l-1,   j-1] = a0p2 + a1p3
                c[l-1,j4+j-1] = a4m6 + a5m7
                c[l-1,j2+j-1] = c2 * a0p2m1p3 - s2 * a4m6m5m7
                c[l-1,j6+j-1] = s2 * a0p2m1p3 + c2 * a4m6m5m7
                c[l-1,j1+j-1] = c1*(a0m2-a5p7)-s1*(a4p6+a1m3)
                c[l-1,j5+j-1] = s1*(a0m2-a5p7)+c1*(a4p6+a1m3)
                c[l-1,j3+j-1] = c3*(a0m2+a5p7)-s3*(a4p6-a1m3)
                c[l-1,j7+j-1] = s3*(a0m2+a5p7)+c3*(a4p6-a1m3)
            jbase=jbase+1
            
        i0 = i0 + iink
        i1 = i1 + iink
        i2 = i2 - iink
        i3 = i3 - iink
        i4 = i4 + iink
        i5 = i5 + iink
        i6 = i6 - iink
        i7 = i7 - iink
        jbase=jbase+7*la

    if (i1 <= i2):
        sin45=np.sqrt(0.5)
        for i in range(1,la+1):
            j=jbase
            for l in range(1,lot+1):
                c[l-1,   j-1]=a[l-1,i0+i-1]+a[l-1,i1+i-1]
                c[l-1,j1+j-1]=sin45*((a[l-1,i0+i-1]-a[l-1,i1+i-1])-(a[l-1,la+i0+i-1]+a[l-1,la+i1+i-1]))
                c[l-1,j2+j-1]=a[l-1,la+i1+i-1]-a[l-1,la+i0+i-1]
                c[l-1,j3+j-1]=-sin45*((a[l-1,i0+i-1]-a[l-1,i1+i-1])+(a[l-1,la+i0+i-1]+a[l-1,la+i1+i-1]))
            jbase=jbase+1
    
    la = la * 4
    
    return c,la


def _ifft3(a,trigs,n,lot):
    
    c = np.zeros((lot,n))
    
    SIN60 = 0.866025403784438

    ib = 2 * (n//3) + 1

    c[:,0] = 0.5 * a[:,0] + a[:,ib-1]
    c[:,1] = 0.5 * a[:,0] - 0.5 * a[:,ib-1] - SIN60 * a[:,ib]
    c[:,2] = 0.5 * a[:,0] - 0.5 * a[:,ib-1] + SIN60 * a[:,ib]

    ia = 3
    ic = ib - 2
    ib = ib + 2

    for j in range(4 , n-7 , 6):
        c1 = trigs[ia-1]
        s1 = trigs[ia  ]
        c2 = trigs[ia+ia-2]
        s2 = trigs[ia+ia-1]

        for l in range(1 , lot+1):
            hbpc = a[l-1,ia-1] - 0.5 * (a[l-1,ib-1] + a[l-1,ic-1])
            hbmc = a[l-1,ia  ] - 0.5 * (a[l-1,ib  ] - a[l-1,ic  ])
            sbmc = SIN60 * (a[l-1,ib-1] - a[l-1,ic-1])
            sbpc = SIN60 * (a[l-1,ib  ] + a[l-1,ic  ])

            c[l-1,j-1] = a[l-1,ia-1] + a[l-1,ib-1] + a[l-1,ic-1]
            c[l-1,j+2] = a[l-1,ia  ] + a[l-1,ib  ] - a[l-1,ic  ]
            c[l-1,j  ] = c1 * (hbpc-sbpc) - s1 * (hbmc+sbmc)
            c[l-1,j+3] = s1 * (hbpc-sbpc) + c1 * (hbmc+sbmc)
            c[l-1,j+1] = c2 * (hbpc+sbpc) - s2 * (hbmc-sbmc)
            c[l-1,j+4] = s2 * (hbpc+sbpc) + c2 * (hbmc-sbmc)
        ia = ia + 2
        ib = ib + 2
        ic = ic - 2


    c[:,n-3] = a[:,ia-1]
    c[:,n-2] =   0.5 * a[:,ia-1] - SIN60 * a[:,ia]
    c[:,n-1] = - 0.5 * a[:,ia-1] - SIN60 * a[:,ia]

    la = 3
    return c,la


def _ifft2(a,trigs,n,lot):
    c = np.zeros((lot,n))

    c[:,0] = 0.5 * a[:,0]
    c[:,1] = c[:,0]

    ia    =   3
    ib    = n-1

    for j in range( 3 , n-4 , 4):
        c1 = trigs[ia-1]
        s1 = trigs[ia  ]
        for l in range(1,lot+1):
            amb = a[l-1,ia-1] - a[l-1,ib-1]
            apb = a[l-1,ia  ] + a[l-1,ib  ]
            c[l-1,j-1] = a[-1,ia-1] + a[l-1,ib-1]
            c[l-1,j+1] = a[-1,ia  ] - a[l-1,ib  ]
            c[l-1,j  ] = c1 * amb - s1 * apb
            c[l-1,j+2] = s1 * amb + c1 * apb
        ia = ia + 2
        ib = ib - 2
    
    c[:,n-2] =  a[:,ia-1]
    c[:,n-1] = -a[:,ia  ]

    la = 2
    return c,la

def _fc2gp(fc,n,lot):
    '''Convert Fourier components to gridpoint space.
    
    Parameters
    ----------
    fc : numpy.ndarray
        Array of Fourier coefficients
    n : int 
        NLON, number of longitudes
    lot : int 
        NLAT*NLEV, number of latitudes times number of longitudes
    '''
    
    a = np.reshape(fc,(lot,n))
    trigs = np.zeros(n) #(n/2-1)*2 + 1 = n-2+1 = n-1 -> len(n)
    
    dell = 4.0 * np.arcsin(1.0)/n
    for k in range(0,n//2):
        angle = k*dell
        trigs[2*k  ] = np.cos(angle)
        trigs[2*k+1] = np.sin(angle)
    
    nf = n//8
    while nf >= 4:
        nf = nf//4
    
    la = 1
    if nf==2:
        a,la = _ifft2(a,trigs,n,lot)
    elif nf==3:
        a,la = _ifft3(a,trigs,n,lot)
    while la < n//8:
        a,la = _ifft4(a,trigs,n,lot,la)
    a = _ifft8(a,n,lot)
    
    return a
        

def readrecord(fbuffer,n,en,ml,mf):
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
    mf : str
        Format of the record marker ('i' or 'l')
        
    Returns
    -------
    array-like, array-like, int
        A tuple containing first the header, then the data contained in the record, and finally the new
        position in the buffer in bytes.
    '''
    if n<len(fbuffer):
        wl,fmt = _getknownwordlength(fbuffer,n,en,ml,mf)
                
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml #Add one word for restatement of header length (for backwards seeking)
        datalength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//wl)
        n+=ml
        data = struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl])
        n+=datalength*wl+ml #additional 4 for restatement of datalength
        return header,data,n
    else:
        raise Exception("Reached end of buffer!!!")
    
def readvariablecode(fbuffer,kcode,en,ml,mf):
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
    mf : str
        Format of the record marker ('i' or 'l')
    
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
        
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml
        datalength = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
        n+=ml
        if header[0]==kcode:
            dataheader = header
            wl, fmt = _getknownwordlength(fbuffer,recordn0,en,ml,mf)
            datalength = int(datalength//wl)
            if not variable:
                variable = np.array(struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            else:
                variable = np.append(variable,struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            n+=datalength*wl+ml
        else: #Fast-forward past this variable without reading it.
            n+=datalength+ml
    
    return dataheader, variable

def _gettimevar(fbuffer):
    '''Extract the time array, as an array of timesteps'''
    
    en = _getEndian(fbuffer)
    ml,mf = _getwordlength(fbuffer,0,en)
    
    kcode = 139 #Use surface temperature to do this
    time = []
    n = 0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml)
    
    variable = None
    
    while n<len(fbuffer):
        
        recordn0 = n
        
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml
        datalength = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
        n+=ml
        if header[0]==kcode:
            time.append(header[6]) #nstep-nstep1 (timesteps since start of run)
        n+=datalength+ml
    
    return time
    
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
    ml,mf = _getwordlength(fbuffer,0,en)
    
    n=0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml,mf)
    
    headers= {'main':mainheader}
    variables = {'main':zsig}
    nlev=mainheader[6]
    variables["sigmah"] = zsig[:nlev]
    variables["time"] = []
    
    while n<len(fbuffer):
        header,field,n = readrecord(fbuffer,n,en,ml,mf)
        kcode = str(header[0])
        if int(kcode)==139:
            variables["time"].append(header[6]) #nstep-nstep1 (timesteps since start of run)
        if kcode not in variables:
            variables[kcode] = np.array(field)
            headers[kcode] = header
            #print("variable %d"%int(kcode))
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
    dim1 = max(header[4],header[5])
    dim2 = min(header[4],header[5])
    if header[1]==1:
        nlevs=nlev
        if len(variable)%(len(variable)//(dim1*dim2*nlevs))!=0:
            nlevs+=1
    else:
        nlevs=1
    ntimes = int(len(variable)//(dim1*dim2*nlevs))
    if nlevs==1:
        if dim2==1:
            newvar = np.reshape(variable,(ntimes,dim1))
        else:
            newvar = np.reshape(variable,(ntimes,dim2,dim1))
    else:
        if dim2==1:
            newvar = np.reshape(variable,(ntimes,nlevs,dim1))
        else:
            newvar = np.reshape(variable,(ntimes,nlevs,dim2,dim1))
            
    return newvar

def readfile(filename):
    '''Extract all variables from a raw plasim output file and refactor them into the right shapes
    
    This routine will only produce what it is in the file; it will not compute derived variables.
    
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
    sigmah = variables['sigmah']
    
    kcodes = list(variables.keys())
    kcodes.remove('sigmah')
    
    #time = _gettimevar(fbuffer) #This means we do one additional sweep through the file
    time = variables["time"]
    kcodes.remove("time")
    
    data = {}
    
    for key in kcodes:
        data[key] = refactorvariable(variables[key],headers[key],nlev=nlevs)
    
    nlat = min(headers['main'][4],headers['main'][5])
    nlon = max(headers['main'][4],headers['main'][5])
    ntru = headers['main'][7]
    ntimes = len(time)
        
    sid,gwd = _gaulat(nlat)
    rlat = np.arcsin(sid)
    lat = rlat*180.0/np.pi
    
    lon = np.arange(nlon)/float(nlon)*360.0
    rlon = lon*np.pi/180.0
    
    sigmab = np.append([0,],sigmah)
    sigma = 0.5*(sigmab[0:-1]+sigmab[1:]) #Mid-layer sigma
    
    data["lev"] = sigma[:]
    data["lat"] = lat[:]
    data["lon"] = lon[:]
    data["time"]= time[:]

    return data

def dataset(filename, variablecodes, mode='grid', zonal=False, substellarlon=0.0, physfilter=False):
    '''Read a raw output file, and construct a dataset.
    
    Parameters
    ----------
    filename : str
        Path to the raw output file
    variablecodes : array-like
        list of variables to include. Can be the integer variable codes from the burn7 postprocessor
        conventions (as either strings or integers), or the short variable name strings 
        (e.g. 'rlut'), or a combination of the two.
    mode : str, optional
        Horizontal output mode. Can be 'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        For grid modes ("grid" and "synchronous"), compute and output zonal means
    substellarlon : float, optional
        If mode='synchronous', the longitude of the substellar point in equatorial coordinates,
        in degrees
    physfilter : bool, optional
        Whether or not a physics filter should be used when transforming spectral variables to
        Fourier or grid domains
        
    Returns
    -------
    dict
        Dictionary of extracted variables
    '''
    
    rawdata = readfile(filename)
    
    
    lat = rawdata["lat"]
    lon = rawdata["lon"]
    lev = rawdata["lev"]
    time = rawdata["time"]
    
    nlat = len(lat)
    nlon = len(lon)
    nlev = len(lev)
    ntime = len(time)
    
    ntru = (nlon-1) // 3
    
    rdataset = {}
    
    for key in variablecodes:
        '''Collect metadata from our built-in list, and extract 
        the variable data if it already exists; if not set a flag
        that we need to derive it.'''
        if type(key)==int:
            try:
                meta = ilibrary[str(key)][:]
            except:
                raise Exception("Unknown variable code requested: %s"%str(key))
            if str(key) in rawdata:
                variable = rawdata[str(key)][:]
                derived=False
            else:
                #print(str(key)+" not in rawdata;",rawdata.keys())
                derived=True
        else:
            if key in ilibrary:
                meta = ilibrary[key][:]
            elif key in slibrary:
                meta = slibrary[key][:]
                kcode = meta[0]
                meta[0] = key
                #print("Reassigning key; key was %s and is now %s"%(key,kcode))
                key = str(kcode) #Now key is always the integer code, and meta[0] is always the name
            else:
                raise Exception("Unknown variable code requested: %s"%key)
            if key in rawdata:
                variable = rawdata[key][:]
                derived=False
            else:
                #print(key+" not in rawdata;",rawdata.keys())
                derived=True
        meta.append(key)
        print(meta,derived,rawdata.keys())
        if not derived:
            #print("Found variable; no need to derive: %s"%meta[0])
            if mode=="grid":
                if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
                    if len(variable.shape)==3: #Include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes*nlevs,variable.shape[2])))
                                   )
                        gridshape = (ntimes,nlevs,nlat,nlon)
                        dims = ["time","lev","lat","lon"]
                    else:
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes,variable.shape[1])))
                                   )
                        gridshape = (ntimes,nlat,nlon)
                        dims = ["time","lat","lon"]
                    gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
                    gridvar = np.reshape(np.transpose(gridvar),gridshape)
                    gridvar = np.roll(gridvar,nlon//2,axis=-1)
                else: #grid variable
                    gridvar = variable
                    if len(variable.shape)==3:
                        dims = ["time","lat","lon"]
                    else:
                        dims = ["time","lev","lat","lon"]
                if meta[0]=="hus":
                    gridvar[gridvar<0] = 0.0
                if zonal:
                    gridvar = np.nanmean(gridvar,axis=-1)
                    dims.remove("lon")
                meta.append(tuple(dims))
                rdataset[meta[0]] = (gridvar,meta)
                
            elif mode=="synchronous":
                if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
                    if len(variable.shape)==3: #Include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes*nlevs,variable.shape[2])))
                                   )
                        gridshape = (ntimes,nlevs,nlat,nlon)
                        dims = ["time","lev","lat","lon"]
                    else:
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes,variable.shape[1])))
                                   )
                        gridshape = (ntimes,nlat,nlon)
                        dims = ["time","lat","lon"]
                    gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
                    gridvar = np.reshape(np.transpose(gridvar),gridshape)
                    gridvar = np.roll(gridvar,nlon//2,axis=-1)
                else: #grid variable
                    gridvar = variable
                    if len(variable.shape)==3:
                        dims = ["time","lat","lon"]
                    else:
                        dims = ["time","lev","lat","lon"]
                if meta[0]=="hus":
                    gridvar[gridvar<0] = 0.0
                if zonal:
                    gridvar = np.nanmean(gridvar,axis=-1)
                    dims.remove("lon")
                lon,lat,tlgridvar = gcmt.eq2tl(gridvar,lon,lat,substellar=substellarlon,
                                               polemethod='interp') #fine bc all vectors are derived
                
                if zonal:
                    tlgridvar = np.nanmean(tlgridvar,axis=-1)
                meta.append(tuple(dims))
                rdataset[meta[0]] = (tlgridvar,meta)
                
            elif mode=="spectral":
                if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
                    specvar = variable
                    if len(variable.shape)==3:
                        dims = ("time","lev","modes")
                    else:
                        dims = ("time","modes")
                else:
                    if len(variable.shape)==4: #Include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        gpvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes*nlevs,nlat,nlon)))
                                   )
                        dims = ("time","lev","modes")
                    else:
                        
                        ntimes = variable.shape[0]
                        gpvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes,nlat,nlon)))
                                   )
                        dims = ("time","modes")
                    spvar = pyfft.gp2sp(gpvar,nlat,nlon,ntru,int(physfilter))
                    specvar = np.transpose(spvar)
                meta.append(dims)
                rdataset[meta[0]] = (specvar,meta)
            elif mode=="fourier":
                if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
                    if len(variable.shape)==3: #Include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes*nlevs,variable.shape[2])))
                                   )
                        fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                        dims = ["time","lev","lat","fourier","complex"]
                    else:
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes,variable.shape[1])))
                                   )
                        fcshape = (ntimes,nlat,nlon//2,2)
                        dims = ["time","lat","fourier","complex"]
                    fcvar = pyfft.sp3fc(spvar,nlat,nlon,ntru,int(physfilter))
                    fouriervar = np.reshape(np.transpose(fcvar),fcshape)
                else: #grid variable
                    if len(variable.shape)==4: #include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        gpvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable/1.4142135623730951,
                                                           (ntimes*nlevs,nlat,nlon)))
                                   )
                        fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                        dims = ["time","lev","lat","fourier","complex"]
                    else:
                        ntimes = variable.shape[0]
                        gpvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable/1.4142135623730951,
                                                           (ntimes,nlat,nlon)))
                                   )
                        fcshape = (ntimes,nlat,nlon//2,2)
                        dims = ["time","lat","fourier","complex"]
                    fcvar = pyfft.gp3fc(gpvar)
                    fouriervar = np.reshape(np.transpose(fcvar),fcshape)
                meta.append(tuple(dims))
                rdataset[meta[0]] = (fouriervar,meta)
                
            elif mode=="syncfourier":
                #First transform into synchronous coordinate space
                if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
                    if len(variable.shape)==3: #Include lev
                        nlevs = variable.shape[1]
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes*nlevs,variable.shape[2])))
                                   )
                        gridshape = (ntimes,nlevs,nlat,nlon)
                        dims = ["time","lev","lat","lon"]
                    else:
                        ntimes = variable.shape[0]
                        spvar = np.asfortranarray(
                                   np.transpose(np.reshape(variable,
                                                           (ntimes,variable.shape[1])))
                                   )
                        gridshape = (ntimes,nlat,nlon)
                        dims = ["time","lat","lon"]
                    gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
                    gridvar = np.reshape(np.transpose(gridvar),gridshape)
                    gridvar = np.roll(gridvar,nlon//2,axis=-1)
                else: #grid variable
                    gridvar = variable
                    if len(variable.shape)==3:
                        dims = ["time","lat","lon"]
                    else:
                        dims = ["time","lev","lat","lon"]
                if meta[0]=="hus":
                    gridvar[gridvar<0] = 0.0
                if zonal:
                    gridvar = np.nanmean(gridvar,axis=-1)
                    dims.remove("lon")
                lon,lat,tlgridvar = gcmt.eq2tl(gridvar,lon,lat,substellar=substellarlon,
                                               polemethod='interp') #fine bc all vectors are derived
                #Next we want to reshape things so that parallel meridians link up to form circles.
                #This will leave us with half the longitudes, and twice the latitudes.
                
                #0 degrees "latitude" fourier coefficients will be the fourier transform along the 
                #substellar polar meridian, while 90 degrees "latitude" will be the transform along
                #the equatorial meridian.

                rottlgridvar = np.zeros(tlgridvar.shape)
                for jlon in range(nlats):
                    rottlgridvar[...,jlon,:nlats] = tlgridvar[...,jlon]
                    rottlgridvar[...,jlon,nlats:] = tlgridvar[...,jlon+nlats]    
                
                #Compute fourier coefficients along our new "longitudes"
                if len(rottlgridvar.shape)==4: #include lev
                    nlevs = rottlgridvar.shape[1]
                    ntimes = rottlgridvar.shape[0]
                    gpvar = np.asfortranarray(
                               np.transpose(np.reshape(rottlgridvar/1.4142135623730951,
                                                       (ntimes*nlevs,nlat,nlon)))
                               )
                    fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                    dims = ["time","lev","lat","fourier","complex"]
                else:
                    ntimes = rottlgridvar.shape[0]
                    gpvar = np.asfortranarray(
                               np.transpose(np.reshape(rottlgridvar/1.4142135623730951,
                                                       (ntimes,nlat,nlon)))
                               )
                    fcshape = (ntimes,nlat,nlon//2,2)
                    dims = ["time","lat","fourier","complex"]
                fcvar = pyfft.gp3fc(gpvar)
                fouriervar = np.reshape(np.transpose(fcvar),fcshape)
                meta.append(tuple(dims))
                rdataset[meta[0]] = (fouriervar,meta)
                
            else:
                raise Exception("Invalid output mode selected")
            
                #### Add options for synchronous fourier and synchronous spectral
                
            print("Collected variable: %s"%meta[0])
        else: #derived=True
            pass            
        
            # Add in derived variables
            #print(nlat,nlon,ntru,nlevs*ntimes)
            #ta = pyfft.sp2gp(np.asfortranarray(np.transpose(np.reshape(data["130"],(ntimes*nlevs,data["130"].shape[-1])))),
                             #nlat,nlon,ntru,int(physfilter))
            #ta = np.reshape(np.transpose(ta),(ntimes,nlevs,nlat,nlon))
            
            #data["ta"] = ta
            
            '''
            burn7 routines we may need to copy:
            * dv2ps (streamfunction etc)
            * spvfc (spectral->fourier)
            * sp2fci (inverse legendre transform)
            * sp2fcd (to get to dpsdy)
            * fc2gp (fourier to gridpoint)
            * Fill humidity holes (where it's -1)
            * Omega_w (to get w wind)
            * sh2rh (relative humidity from specific humidity)
            * Extrap (to get sea-level pressure)
            * Speed (for windspeed)
            * Remember to add 142 and 143 for precip 
            * Add 178 and 179 for net_top 
            * add 176 and 177 for net_bot 
            * Add 218*L_times_rhoH2O and 176 and 177 and 146 and 147 for net_heat
            * add 182 - 160 + 142 + 143 for net_water 
            * sw_atm -> 178 - 176
            * lw_atm -> 179 - 177
            * net_atm -> 178 + 179 - 176 - 177
            * surf_runoff -> 182 - 221 + 142 + 143
            * fresh_water -> 142 + 143 + 182
            * gp2fc_uv
            * fc2sp_uv
            * sp2fc_uv
            * fc2gp_uv
            * gp2fc 
            * scaluv
            * fc2sp 
            '''
                
        
          
    rdataset["lat"] = (lat,["lat","latitude","deg"])
    rdataset["lon"] = (lon,["lon","longitude","deg"])
    rdataset["lev"] = (lev,["lev","sigma_coordinate","nondimensional"])
    rdataset["time"] = (time,["time","timestep_of_year","timesteps"])      
    
    return rdataset

                
def netcdf(rdataset,filename="most_output.nc"):
    '''Write a dataset to a netCDF file.
    
    Parameters
    ----------
    dataset : dict
        A dictionary of outputs as generated from :py:func`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written.
        
    Returns
    -------
    object
        A netCDF object corresponding to the file that has been written.
    '''
    
    ncd = nc.Dataset(filename, "w", format="NETCDF4")
    
    latitude  = rdataset["lat" ]
    longitude = rdataset["lon" ]
    level     = rdataset["lev" ]
    timestamp = rdataset["time"]
    
    nlats  = len( latitude[0])
    nlons  = len(longitude[0])
    nlevs  = len(    level[0])
    ntimes = len(timestamp[0])
    ntru = (nlons-1)//3
    nmodes = (ntru+1)*(ntru+2)
    
    lat = ncd.createDimension("lat",   nlats)
    lon = ncd.createDimension("lon",   nlons)
    lev = ncd.createDimension("lev",   nlevs)
    ttime = ncd.createDimension("time",ntimes)
    cmplx = ncd.createDimension("complex",2)
    fourier = ncd.createDimension("fourier",nlons//2)
    sphmods = ncd.createDimension("modes",nmodes)
    
    latitudes   = ncd.createVariable("lat", "f4",("lat", ),zlib=True,least_significant_digit=6)
    longitudes  = ncd.createVariable("lon", "f4",("lon", ),zlib=True,least_significant_digit=6)
    levels      = ncd.createVariable("lev", "f4",("lev", ),zlib=True,least_significant_digit=6)
    times       = ncd.createVariable("time","f4",("time",),zlib=True,least_significant_digit=6)
    complexn    = ncd.createVariable("complex","f4",("complex",),zlib=True,least_significant_digit=6)
    fourierc    = ncd.createVariable("fharmonic","f4",("fourier",),zlib=True,least_significant_digit=6)
    spharmonics = ncd.createVariable("modes","f4",("modes",),zlib=True,least_significant_digit=6)
    
    ncd.set_auto_mask(False)
    latitudes.set_auto_mask(False)  
    longitudes.set_auto_mask(False)   
    levels.set_auto_mask(False)
    times.set_auto_mask(False)
    complexn.set_auto_mask(False)
    fourierc.set_auto_mask(False)
    spharmonics.set_auto_mask(False)
    
    latitudes.units  =  latitude[1][2]
    longitudes.units = longitude[1][2]
    levels.units     =     level[1][2]
    times.units      = timestamp[1][2]
    complexn.units   = "n/a"
    fourierc.units   = "n/a"
    spharmonics.units= "n/a"
    
    latitudes[:]   =  latitude[0].astype("float32")
    longitudes[:]  = longitude[0].astype("float32")
    levels[:]      =     level[0].astype("float32")
    times[:]       = timestamp[0].astype("float32")
    complexn[0] = np.float32(1.0)
    complexn[1] = np.float32(1.0j)
    fourierc[:] = np.arange(nlons//2,dtype='float32')
    sphmods[:] = np.arange(nmodes,dtype='float32')
    
    longitudes.axis = 'X'
    fourierc.axis   = 'X'
    sphmods.axis    = 'X'
    latitudes.axis  = 'Y'
    levels.axis     = 'Z'
    
    latitudes.standard_name = latitude[1][1]
    longitudes.standard_name= longitude[1][1]
    levels.standard_name    = levels[1][1]
    complexn.standard_name  = "complex_plane"
    fourierc.standard_name  = "fourier_coefficients"
    sphmods.standard_name   = "spherical_real_modes"
    times.standard_name     = timestamp[1][1]
    
    latitudes.long_name = latitude[1][1]
    longitudes.long_name= longitude[1][1]
    levels.long_name    = "sigma at layer midpoints"
    complexn.long_name  = "complex coefficients"
    fourierc.long_name  = "Fourier coefficients"
    sphmods.long_name   = "Spherical harmonic real global modes"
    times.long_name     = timestamp[1][1]
    
    levels.positive = "down"
    
    keyvars = rdataset.keys()
    keyvars.remove("time")
    keyvars.remove("lat" )
    keyvars.remove("lon" )
    keyvars.remove("lev" )
    
    for key in keyvars:
        datavar,meta = rdataset[key]
        dims = meta[4]
        
        variable = ncd.createVariable(key,"f4",dims,zlib=True,least_significant_digit=6)
        variable.set_auto_mask(False)
        variable.units = meta[2]
        variable[:] = datavar[:]
        variable.standard_name = meta[1]
        variable.long_name = meta[1]
        variable.units = meta[2]
        variable.code = meta[3]
        if "fourier" not in dims and "modes" not in dims:
            variable.grid_type = "gaussian"
        print("Packing %s in %s\t....... %d timestamps"%(key,filename,ntimes))
        
    ncd.sync()
    return ncd
    
#Constants

MARS_GRAV   = 3.728
MARS_RADIUS = 3400000.0
MARS_RD     = 189.0


def postprocess(rawfile,outfile,namelist=None,variables=ilibrary.keys(),dformat="netcdf",mode='grid',
                zonal=False, substellarlon=0.0, physfilter=False.timeaverage=True,stdev=False
                times=12,interpolatetimes=True,radius=1.0,gravity=9.80665,gascon=287.0):
    '''Convert a raw output file into a postprocessed formatted file.
    
    Parameters
    ----------
    rawfile : str
        Path to the raw output file
    outfile : str
        Path to the destination output file
    namelist : str, optional
        Path to a burn7 postprocessor namelist file. If not given, then `variables` must be set. 
    variables : list or dict, optional
        If a list is given, a list of either variable keycodes (integers or strings), or the abbreviated
        variable name (e.g. 'ts' for surface temperature). If a dict is given, each item in the dictionary
        should have the keycode or variable name as the key, and the desired horizontal mode and additional
        options for that variable as a sub-dict. Each member of the subdict should be passable as **kwargs 
        to :py:func`advancedDataset() <exoplasim.pyburn.advancedDataset>`. If None, then `namelist` must be set.
    dformat : str, optional
        Output format. At present, only "netcdf" and "numpy" (numpy.savez) are supported.
    mode : str, optional
        Horizontal output mode, if modes are not specified for individual variables. Options are 
        'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        Whether zonal means should be computed for applicable variables.
    substellarlon : float, optional
        Longitude of the substellar point. Only relevant if a synchronous coordinate output mode is chosen.
    physfilter : bool, optional
        Whether or not a physics filter should be used in spectral transforms.
    times : int or array-like, optional
        Either the number of timestamps by which to divide the output, or a list of times given as a fraction
        of the output file duration (which enables e.g. a higher frequency of outputs during periapse of an
        eccentric orbit, when insolation is changing more rapidly).
    timeaverage : bool, optional
        Whether or not timestamps in the output file should be averaged to produce the requested number of 
        output timestamps. Timestamps for averaged outputs will correspond to the middle of the averaged time period.
    stdev : bool, optional
        Whether or not standard deviations should be computed. If timeaverage is True, this will be the 
        standard deviation over the averaged time period; if False, then it will be the standard deviation
        over the whole duration of the output file
    interpolatetimes : bool, optional
        If true, then if the times requested don't correspond to existing timestamps, outputs will be
        linearly interpolated to those times. If false, then nearest-neighbor interpolation will be used.
    radius : float, optional
        Planet radius in Earth radii
    gravity : float, optional
        Surface gravity in m/s^2.
    gascon : float, optional
        Specific gas constant for dry gas (R$_d$) in J/kg/K.  
        
    
    '''
    
    if namelist is None:
        if variables is None:
            raise Exception("Must specify either a burn7 namelist or provide a list or dict of variables")
        if type(variables)==tuple or type(variables)==list:
            data = dataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,gascon=gascon, 
                           zonal=zonal, substellarlon=substellarlon, physfilter=physfilter)
        elif type(variables)==dict: #for advancedDataset
            pass
            #data = advancedDataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,
                                   #gascon=gascon, zonal=zonal, substellarlon=substellarlon, 
                                   #physfilter=physfilter)
        
    else:
        #Scrape namelist
        variables=None
        with open(namelist,"r") as fname:
            fnamelist = fname.read().split('\n')
        for line in fnamelist:
            parts = line.split('=')
            key = parts[0]
            arg = parts[1]
            if key=="code" or key=="CODE":
                variables = arg.split(',')
            elif key=="htype" or key=="HYTPE":
                if arg=="g" or arg=="G":
                    mode="grid"
                elif arg=="z" or arg=="Z":
                    mode="grid"
                    zonal=True
                elif arg=="s" or arg=="S":
                    mode="spectral"
                elif arg=="f" or arg=="F":
                    mode="fourier"
            elif key=="mean" or key=="MEAN":
                if int(arg)==0:
                    timeaverage=False
                    stdev = False
                elif int(arg)==1:
                    timeaverage=True
                    stdev = False
                elif int(arg)==2:
                    stdev = True
                elif int(arg)==3:
                    timeaverage=True
                    stdev = True
            elif key=="radius" or key=="RADIUS":
                radius = float(arg)/6371220.0      #Radii in the namelist are given in metres
            elif key=="gravity" or key=="GRAVITY":
                gravity = float(arg)
            elif (key=="mars" or key=="MARS") and int(arg)==1:
                gravity = MARS_GRAV
                radius  = MARS_RADIUS/6371220.0    #We want to start off with radii in Earth radii
                gascon  = MARS_RD      #This is called RD in burn7, not gascon
        data = dataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,gascon=gascon, 
                       zonal=zonal, substellarlon=substellarlon, physfilter=physfilter)
        
    # Compute time averages, binning, stdev, etc
    
    dtimes = data["time"][0]
    ntimes = len(dtimes)
    
    if type(times)==int: #A number of time outputs was specified
        if ntimes==times: #The number of outputs exactly equals the number provided
            pass #We don't need to do anything
        
        else:
            if timeaverage:
               indices = np.linspace(0,ntimes,times+1,True).astype(int)
               counts = np.diff(indices)
               varkeys = data.keys()
               varkeys.remove("time")
               varkeys.remove("lat")
               varkeys.remove("lon")
               varkeys.remove("lev")
               newtimes = np.add.reduceat(dtimes,indices[:-1]) / counts
               for var in varkeys:
                   odata = data[var][0][:]
                   newshape = list(odata.shape)
                   newshape[0] = times
                   newshape = tuple(newshape)
                   data[var][0] = np.add.reduceat(odata,indices[:-1],axis=0) / np.resize(counts,newshape)
               data["time"][0] = newtimes
            else:
               if interpolatetimes:
                   interpolation = "linear"
               else:
                   interpolation = "nearest"
               newtimes = np.linspace(0.,1.,num=times)*dtimes[-1] #We can't extrapolate; only interpolate
               varkeys = data.keys()
               varkeys.remove("time")
               varkeys.remove("lat")
               varkeys.remove("lon")
               varkeys.remove("lev")
               for var in varkeys:
                   odata = data[var][0][:]
                   interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind=interpolation)
                   data[var][0] = interpfunc(newtimes)
               if not interpolatetimes:
                   interpfunc = scipy.interpolate.interp1d(dtimes,dtimes,kind="nearest")
                   data["time"][0] = interpfunc(newtimes)
               else:
                   data["time"][0] = newtimes
                   
    else: #A list of times was specified
        
        if timeaverage: #times values are assumed to be bin edges
            if not interpolatetimes: #We will always round down to the nearest neighbor
                indices = np.digitize(np.array(times)*dtimes[-1],dtimes)-1
                counts = np.diff(indices)
                varkeys = data.keys()
                varkeys.remove("time")
                varkeys.remove("lat")
                varkeys.remove("lon")
                varkeys.remove("lev")
                newtimes = np.add.reduceat(dtimes,indices[:-1]) / counts
                for var in varkeys:
                    odata = data[var][0][:]
                    newshape = list(odata.shape)
                    newshape[0] = times
                    newshape = tuple(newshape)
                    data[var][0] = np.add.reduceat(odata,indices[:-1],axis=0) / np.resize(counts,newshape)
                data["time"][0] = newtimes
            else: #First we'll interpolate to high time resolution, then compute average via binning
                ttimes = np.linspace(dtimes[0],dtimes[-1],num=10*ntimes)
                indices = np.digitize(np.array(times)*dtimes[-1],ttimes)-1
                counts = np.diff(indices)
                varkeys = data.keys()
                varkeys.remove("time")
                varkeys.remove("lat")
                varkeys.remove("lon")
                varkeys.remove("lev")
                for var in varkeys:
                    odata = data[var][0][:]
                    interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind='linear')
                    tempdata = interpfunc(ttimes)
                    newshape = list(odata.shape)
                    newshape[0] = len(times)
                    newshape = tuple(newshape)
                    data[var][0] = np.add.reduceat(tempdata,indices[:-1],axis=0) / np.resize(counts,newshape)
                newtimes = np.array(times)
                newtimes = 0.5*(newtimes[:-1]+newtimes[1:])*dtimes[-1]
                data["time"][0] = newtimes
                    
        else:
            newtimes = np.array(times)*dtimes[-1] #Convert to timestamps
            if interpolatetimes:
                interpolation = "linear"
            else:
                interpolation = "nearest"
            varkeys = data.keys()
            varkeys.remove("time")
            varkeys.remove("lat")
            varkeys.remove("lon")
            varkeys.remove("lev")
            for var in varkeys:
                odata = data[var][0][:]
                interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind=interpolation)
                data[var][0] = interpfunc(newtimes)
            if not interpolatetimes:
                interpfunc = scipy.interpolate.interp1d(dtimes,dtimes,kind="nearest")
                data["time"][0] = interpfunc(newtimes)
            else:
                data["time"][0] = newtimes
        
    #### Need to add standard deviation calculation; should go in the timeaverage parts of the computation
    #### above, or for when timeaverage is False, it can go here (then we compute the std dev of the whole
    #### timeseries)
    
        
    
        
    
    
    # Write to output
    
    if dtype=="netcdf": #outfile.nc
        netcdf(data,filename=outfile)
        
    elif dtype=="numpy": #outfile.npz
        npsavez(data,filename=outfile)
        
        
        
    
    
