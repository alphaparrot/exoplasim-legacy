================================
Postprocessing ExoPlaSim Outputs
================================

The Basics: Formats, Variables, and Math
----------------------------------------

As of ExoPlaSim 3.0.0, postprocessing can be done using the :py:mod:`exoplasim.pyburn <exoplasim.pyburn>`
module. This module exposes an API for setting the variables to be included in postprocessed output,
the horizontal mode in which to present them, and any additional math that should be performed, including
coordinate transformations, time-averaging, and standard deviations. ``pyburn`` also supports a large
range of output formats: netCDF, HDF5, NumPy's compressed ``.npz`` archives (the default), and plain-text
comma-separated value (CSV) files. The latter can be compressed individually with the ``gzip`` format,
tarballed, or tarballed and compressed (in the latter case with ``gzip``, ``lzma``, or ``bzip2`` 
compression types). Producing netCDF files requires that the netCDF4 python library be present (you 
can install it with ``pip install netCDF4`` or at ExoPlaSim's install-time with ``pip install exoplasim[netCDF4]``). Similarly, producing HDF5 files requires the presence of the ``h5py`` Python
library, which can be installed via ``pip install h5py`` or, at ExoPlaSim's install time, with
``pip install exoplasim[HDF5]``. Support for both netCDF and HDF5 can be guaranteed at install-time
by combining them::

    pip install exoplasim[netCDF4,HDF5]

Format
******
    
The choice of output format can be specified either when the postprocessor is called (if being used
manually), or as an argument to a :py:class:`Model <exoplasim.Model>` object, by simply providing the
file extension:

    +-----------------+----------------------+
    | Format          | Supported Extensions |
    +=================+======================+
    | NumPy (default) | .npz                 |
    |                 +----------------------+
    |                 | .npy                 |
    +-----------------+----------------------+
    | netCDF          | .nc                  |
    +-----------------+----------------------+
    | HDF5            | .h5                  |
    |                 +----------------------+
    |                 | .he5                 |
    |                 +----------------------+
    |                 | .hdf5                |
    +-----------------+----------------------+
    | Compressed CSV  | .gz                  |
    |                 +----------------------+
    |                 | .tar.gz              |
    |                 +----------------------+
    |                 | .tar.xz              |
    |                 +----------------------+
    |                 | .tar.bz2             |
    +-----------------+----------------------+
    | Uncompressed CSV| .csv                 |
    |                 +----------------------+
    |                 | .txt                 |
    |                 +----------------------+
    |                 | .tar                 |
    +-----------------+----------------------+
    
Because the NumPy archive format does not support additional metadata arrays, metadata is stored
separately in a file using the ``_metadata.npz`` suffix. This file is typically a few tens of kiB. 

CSV-type files will only contain 2D variable information, so the first N-1 dimensions will be flattened. 
The original variable shape is included in the file header (prepended with a # character) as the first 
items in a comma-separated list, with the first non-dimension item given as the '|||' placeholder. On 
reading variables from these files, they should be reshaped according to these dimensions. This is true 
even in tarballs (which contain CSV files). If read in by :py:func:`gcmt.load() <exoplasim.gcmt.load>`, 
this reshaping will be done automatically.

Note that when using the :py:func:`pyburn.postprocess() <exoplasim.pyburn.postprocess>` function
directly, **a single file must be specified as the output file.** This is true even for formats
that produce a large number of files that don't get bound up together, such as ``.gz`` and ``.csv``,
which produce a folder containing one file per variable. The file you specify should have the pattern
``<subdirectory>.<extension>``. This file will not actually be created, but it will be parsed to
determine the desired output format. So, for example, to create an archive consisting of a folder full
of CSV files for the raw output file ``MOST.00127``, one would use ``MOST.00127.csv``. The surface
temperature variable, ``ts``, would then be found in ``MOST.00127/MOST.00127_ts.csv``. 
**This same combined-format fictional filestring can be passed to **
:py:func:`gcmt.load() <exoplasim.gcmt.load>`**.** The object returned by that function will access the
data in the archive just as if it were a bound archive, such as a tarball, netCDF file, or HDF5 file.
    
A T21 model output with 10 vertical levels, 12 output times, all supported variables in grid 
mode,and no standard deviation computation will have the following sizes for each format:
    
    +-----------------+-----------+
    | Format          | Size      |
    +=================+===========+
    | netCDF          | 12.8 MiB  |
    +-----------------+-----------+
    | HDF5            | 17.2 MiB  |
    +-----------------+-----------+
    | NumPy (default) | 19.3 MiB  |
    +-----------------+-----------+
    | tar.xz          | 33.6 MiB  |
    +-----------------+-----------+
    | tar.bz2         | 36.8 MiB  |
    +-----------------+-----------+
    | gzipped         | 45.9 MiB  |
    +-----------------+-----------+
    | uncompressed    | 160.2 MiB |
    +-----------------+-----------+
            
Variables
*********

Output variables can be chosen in multiple ways. Either a ``burn7``-style namelist can be provided,
containing a list of numeric variables codes (listed below), or a list can be passed directly, containing
a list of numeri codes, a list of strings of numeric codes, or a list of string variable keys, as 
indicated in the leftmost-column of the table below.

Variable lists can be specified once for all outputs of a given type ('regular', 'snapshot', or 
'highcadence'), with :py:func:`Model.cfgpostprocessor() <exoplasim.Model.cfgpostprocessor>`, or
for each model year with :py:func:`Model.postprocess() <exoplasim.Model.postprocess>`, or manually
outside of the ExoPlaSim :py:class:`Model <exoplasim.Model>` object, with 
:py:func:`pyburn.postprocess <exoplasim.pyburn.postprocess>`.

Optionally, as advanced usage, a dictionary can be passed, with one member per variable (using the same
identification rules given above), and :py:func:`pyburn.dataset() <exoplasim.pyburn.dataset>`
keyword arguments specified for each variable. For example, to create an output file with two variables,
surface temperature and streamfunction, both on a horizontal grid, and the streamfunction 
zonally-averaged and passed through physics filters::

    {"ts":{"mode":"grid","zonal":False},
     "stf":{"mode":"grid","zonal":True,"physfilter":True}}

This can be specified in one of 3 ways. Either it can be set for all outputs of a given type
('regular', 'snapshot', or 'highcadence') as a Model property:

>>> myModel.cfgpostprocessor(ftype="regular",extension=".nc",
>>>                          variables={"ts":{"mode":"grid","zonal":False},
>>>                                     "stf":{"mode":"grid","zonal":True,"physfilter":True}})
                                    
Or it can be set each time :py:func:`Model.postprocess() <exoplasim.Model.postprocess>` is called:

>>> myModel.postprocess("MOST.00127",
>>>                     {"ts":{"mode":"grid","zonal":False},
>>>                      "stf":{"mode":"grid","zonal":True,"physfilter":True}},
>>>                     log="burnlog.00127",crashifbroken=True)
                        
Or, finally, it can be specified directly to 
:py:func:`pyburn.postprocess() <exoplasim.pyburn.postprocess>`:

>>> pyburn.postprocess("MOST.00127","MOST.00127.nc",logfile="burnlog.00127",
>>>                    variables={"ts":{"mode":"grid","zonal":False},
>>>                               "stf":{"mode":"grid","zonal":True,"physfilter":True}})
                                  
Postprocessing Math
*******************

``pyburn`` provides the ability to perform various mathematical operations on the data as part of
the postprocessing step. 

Multiple horizontal modes are available (specified with the ``mode`` keyword), including a 
Gaussian-spaced latitude-longitude grid (``"grid"``), spherical harmonics (``"spectral"``), 
Fourier coefficients for each latitude (``"fourier"``), a latitude-longitude grid rotated such that the 
"North" pole is at the substellar point of a sychronously-rotating planet, and the "equator" is the 
terminator (``"synchronous"``), and Fourier coefficients computed along lines of constant longitude 
(including the mirror component on the opposite hemisphere) in that rotated coordinate system 
(``"syncfourier"``). Additionally, for output modes with discrete latitudes, variables can be 
zonally-averaged (``zonal=True``).

ExoPlaSim performs some time-averaging on the fly (for "regular"-type outputs) to avoid overloading 
I/O buffers and creating enormous raw output files, but the number of output times is still often
going to be more than you prefer for the postprocessed output data. The default configuration,
for example, produces 72 output timestamps per year. ``pyburn`` can perform time-averaging to reduce
this to e.g. monthly output, via the ``times`` keyword and the ``timeaveraging`` keyword. The former
specifies either the number of output times or the specific output times requested (as decimal fractions
of a model output's timeseries), while the latter is a boolean True/False flag. If specific output times
are requested or the number of requested outputs doesn't divide cleanly into the number of timestamps
in the raw output, ``pyburn`` can interpolate between timestamps using linear interpolation. No 
extrapolation is performed, so you cannot request a time between e.g. the last output of the previous
year and the first output of the current year. Whether or not linear interpolation is used or 
"nearest-neighbor" interpolation (which simply selects the timestamp closest to the target time) can 
be set with the ``interpolatetimes`` keyword--if ``True``--linear interpolation will be used when 
necessary. The minimum number of timestamps in the output file is 1; this corresponds to computing an
annual average. 

Finally, ``pyburn`` brings the ability to compute the standard deviations of ExoPlaSim variables.
Enabling this with ``stdev=True`` will compute the standard deviation in one of two ways, depending
on whether time-averaging is being used. If time averages are being computed, then a standard deviation
will be computed **alongside** each average, and the each standard deviation variable (denoted with the
``_std`` suffix in the variable name, e.g. ``ts_std`` for the standard deviation of surface temperature)
will have the same number of timestamps as the time-averages. If time-averages are **not** being 
computed, then the standard deviation of the entire file's timeseries will be computed, and there will
be one timestamp per standard deviation variable.

Reading Postprocessed Files
---------------------------

While postprocessed files are portable and can be read however you like, ExoPlaSim also provides a
native, format-agnostic way to access them via the :py:func:`gcmt.load() <exoplasim.gcmt.load>`
function. This takes the archive filename as its argument, and returns an object analogous to an 
open netCDF file object. It has two members of interest to the user: ``variables`` and ``metadata``.
Both are compatible with all dictionary methods, and individual variables' data can be extracted by
using the variable name as the dictionary key. For example:

>>> import exoplasim.gcmt as gcmt
>>> myData = gcmt.load("MOST.0127.tar.gz")
>>> surfacetemperature = myData.variables['ts']
>>> surftemp_metadata  = myData.metadata['ts'] 
    
Note that for CSV-type formats, like the tarball given above, the file is left compressed (except
during the initial read), and the whole dataset is `not` loaded into memory. Dimension arrays,
such as latitude, longitude, etc, are loaded, as is all metadata. By default, however, only one
data array will be loaded into memory. This can be expanded with the ``csvbuffersize`` keyword,
which takes the number of variables to permit to hold in the memory buffer. This buffer uses a
first-in, first-out approach, so if a new variable is requested and the buffer is full, the loaded
variable which was used the least recently will be purged from memory.

Postprocessor Variable Codes
----------------------------

Note that in addition to the variable codes listed below, if ``pyburn`` is used with ``stdev=True``,
there will also be variables that correspond to those listed below, with the ``_std`` suffix. If
time-averaging was performed during postprocessing, the standard deviation will be the standard deviation
for each averaged time period, and there will be the same number of timestamps for the ``_std`` variables
as for their nominal data counterparts. If time-averaging was not used, then each standard deviation
variable will have only one timestamp, corresponding to the standard deviation throughout the entire
timeseries present in the file.

+----------+-------+----------------------------------------+---------------------------+---------------+
| Variable |  Code |  Description                           | Units                     |  Notes        |
+==========+=======+========================================+===========================+===============+
|   mld    |  110  |  mixed layer depth                     |  m                        |               |  
+----------+-------+----------------------------------------+---------------------------+---------------+
|   sg     |  129  |  surface geopotential                  |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ta     |  130  |  air temperature                       |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ua     |  131  |  eastward wind                         |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   va     |  132  |  northward wind                        |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hus    |  133  |  specific humidity                     |  kg/kg                    |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ps     |  134  |  surface air pressure                  |  hPa                      |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   wap    |  135  |  vertical air velocity                 |  Pa s-1                   |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   wa     |  137  |  upward wind                           |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   zeta   |  138  |  atm relative vorticity                |  s\ :sup:`-1`             |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ts     |  139  |  surface temperature                   |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   mrso   |  140  |  lwe of soil moisture content          |  m                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   snd    |  141  |  surface snow thickness                |  m                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   prl    |  142  |  lwe of large scale precipitation      |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   prc    |  143  |  convective precipitation rate         |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   prsn   |  144  |  lwe of snowfall amount                |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   bld    |  145  |  dissipation in boundary layer         |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hfss   |  146  |  surface sensible heat flux            |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hfls   |  147  |  surface latent heat flux              |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   stf    |  148  |  streamfunction                        |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   psi    |  149  |  velocity potential                    |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   psl    |  151  |  air pressure at sea level             |  hPa                      |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   pl     |  152  |  log surface pressure                  |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   d      |  155  |  divergence of wind                    |  s\ :sup:`-1`             |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   zg     |  156  |  geopotential height                   |  m                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hur    |  157  |  relative humidity                     |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   tps    |  158  |  tendency of surface air pressure      |  Pa s-1                   |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   u3     |  159  |  u*                                    |  m\ :sup:`3` s\ :sup:`-3` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   mrro   |  160  |  surface runoff                        |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   clw    |  161  |  liquid water content                  |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   cl     |  162  |  cloud area fraction in layer          |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   tcc    |  163  |  total cloud cover                     |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   clt    |  164  |  cloud area fraction                   |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   uas    |  165  |  eastward wind 10m                     |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   vas    |  166  |  northward wind 10m                    |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tas    |  167  |  air temperature 2m                    |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   td2m   |  168  |  dew point temperature 2m              |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tsa    |  169  |  surface temperature accumulated       |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tsod   |  170  |  deep soil temperature                 |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dsw    |  171  |  deep soil wetness                     |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   lsm    |  172  |  land binary mask                      |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   z0     |  173  |  surface roughness length              |  m                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   alb    |  174  |  surface albedo                        |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   as     |  175  |  surface albedo                        |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   rss    |  176  |  surface net shortwave flux            |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   rls    |  177  |  surface net longwave flux             |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   rst    |  178  |  toa net shortwave flux                |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   rlut   |  179  |  toa net longwave flux                 |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+ 
|   tauu   |  180  |  surface eastward stress               |  Pa                       |               |
+----------+-------+----------------------------------------+---------------------------+---------------+ 
|   tauv   |  181  |  surface northward stress              |  Pa                       |               |
+----------+-------+----------------------------------------+---------------------------+---------------+ 
|   evap   |  182  |  lwe of water evaporation              |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+ 
|   tso    |  183  |  climate deep soil temperature         |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   wsoi   |  184  |  climate deep soil wetness             |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   vegc   |  199  |  vegetation cover                      |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   rsut   |  203  |  toa outgoing shortwave flux           |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ssru   |  204  |  surface solar radiation upward        |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   stru   |  205  |  surface thermal radiation upward      |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tso2   |  207  |  soil temperature level 2              |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tso3   |  208  |  soil temperature level 3              |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   tso4   |  209  |  soil temperature level 4              |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   sic    |  210  |  sea ice cover                         |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   sit    |  211  |  sea ice thickness                     |  m                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   vegf   |  212  |  forest cover                          |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   snm    |  218  |  snow melt                             |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   sndc   |  221  |  snow depth change                     |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   prw    |  230  |  atmosphere water vapor content        |  kg m\ :sup:`-2`          |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   glac   |  232  |  glacier cover                         |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+       
|   tsn    |  238  |  snow temperature                      |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   spd    |  259  |  wind speed                            |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   pr     |  260  |  total precipitation                   |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   ntr    |  261  |  net top radiation                     |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   nbr    |  262  |  net bottom radiation                  |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hfns   |  263  |  surface downward heat flux            |  W m\ :sup:`-2`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   wfn    |  264  |  net water flux                        |  m s\ :sup:`-1`           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   lwth   |  266  |  local weathering                      |  W earth                  |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   grnz   |  267  |  ground geopotential                   |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   icez   |  301  |  glacier geopotential                  |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   netz   |  302  |  net geopotential                      |  m\ :sup:`2` s\ :sup:`-2` |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dpdx   |  273  |  d(ps)/dx                              |  Pa m\ :sup:`-1`          |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dpdy   |  274  |  d(ps)/dy                              |  Pa m\ :sup:`-1`          |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   hlpr   |  277  |  half level pressure                   |  Pa                       |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   flpr   |  278  |  full level pressure                   |  Pa                       |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   thetah |  279  |  half level potential temperature      |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   theta  |  280  |  full level potential temperature      |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   czen   |  318  |  cosine solar zenith angle             |  nondimensional           |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   wthpr  |  319  |  weatherable precipitation             |  mm day\ :sup:`-1`        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   mint   |  320  |  minimum temperature                   |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   maxt   |  321  |  maximum temperature                   |  K                        |               |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   cape   |  322  |  convective available potential energy |  J kg\ :sup:`-1`          |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   lnb    |  323  |  level of neutral buoyancy             |  hPa                      |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   sdef   |  324  |  troposphere entropy deficit           |  nondimensional           |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+      
|   absz   |  325  |  sigma-0.85 abs vorticity              |  s\ :sup:`-1`             |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   umax   |  326  |  maximum potential intensity           |  m s\ :sup:`-1`           |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   vent   |  327  |  ventilation index                     |  nondimensional           |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+      
|   vrumax |  328  |  ventilation-reduced maximum wind      |  m s\ :sup:`-1`           |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   gpi    |  329  |  genesis potential index               |  nondimensional           |  Storm Clim.  |
+----------+-------+----------------------------------------+---------------------------+---------------+   
|   dfu    |  404  |  shortwave up                          |  W m\ :sup:`-2`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dfd    |  405  |  shortwave down                        |  W m\ :sup:`-2`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dftu   |  406  |  longwave up                           |  W m\ :sup:`-2`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dftd   |  407  |  longwave down                         |  W m\ :sup:`-2`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dtdt   |  408  |  radiative heating rate                |  K s\ :sup:`-1`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
|   dfdz   |  409  |  flux convergence                      |  W m\ :sup:`-3`           | Snapshot Only |
+----------+-------+----------------------------------------+---------------------------+---------------+
 
 
Burn7 Postprocessor Options
---------------------------

The C++ ``burn7`` postprocessor is now deprecated. It is still included in the model as a legacy option,
and as a backup in case of problems with the ``pyburn`` postprocessor, but as it depends on legacy,
unsupported netCDF-C and netCDF-C++ libraries, compilation can be tricky. For the most part, ``pyburn``
provides all the same features and more, with smaller output files, more output formats, and more
postprocessor options. Documentation for the ``burn7`` postprocessor is however included below for
legacy users who wish to continue using ``burn7``. To use ``burn7`` instead of ``pyburn``, pass
``burn7=True`` when instantiating a :py:class:`Model <exoplasim.Model>`.

Building and Compiling
**********************

To build ``burn7`` for the very first time, in the exoplasim/postprocessor folder, run:::

    ./build_init.sh
    
**You must have a netCDF C++ library already available, which contains netcdf.h.** A
patched NetCDF library will be built with the features burn7 needs, and it will need
your netcdf.h.

To recompile the postprocessor without rebuilding the patched NetCDF library, in the 
exoplasim/postprocessor folder, run:::
    
    ./compile.sh
    
Manual Usage
************

To postprocess raw ExoPlaSim output into NetCDF files, you will need to provide ``burn7``
with a namelist, direct its standard output to a log file, specify the raw output file,
and specify an output file (which does not need to exist yet):::

    ./burn7.x -n<example.nl>burnout [raw_output_file] [output_file]
    
If instead of ``-n`` you use ``-g``, ``.srv`` SERVICE files will be written, and a GraDS
control file will be create--e.g. if ``output_file`` is ``mymodel``, the output will be 
``mymodel.srv`` and ``mymodel.ctl``. If you wish to use NetCDF output (the default), 
you will need to give ``mymodel.nc`` as ``output_file`` instead of simply ``mymodel``.`

The namelist file can take the following parameters:

Parameters
##########
    code : list of integers
        A comma-separated list of integers, corresponding to the output codes listed in the table below.
    vtype : {S, P}
        Whether output vertical structure should be by sigma level (as is native for PlaSim), or converted to pressure levels. Strongly recommend keeping this to S, for sigma levels.
    htype : {S, F, Z, G}
        Whether horizontal representation should be as spherical harmonics (S), Fourier coefficients (F), zonal averagez (Z), or Gauss-spaced latitude-longitude grid (G). Strongly recommend keeping this to G, unless you know what you're doing.
    mean : {0, 1, 2, 3}
        What kind of averaging to do. 3 is the default, and 0 is no time averaging. 1 is monthly means, and 2 is monthly standard deviations. 3 is a combination of both. Mean!=0 can only be used with htype=G.
    netcdf : {0, 1}
        If 1, write NetCDF output file; if 0, write ``.srv`` SERVICE output. 
    mars : {0, 1} (optional)
        If 1, use Mars defaults for computing derived variables.
    gravity : float (optional)
        Surface gravity in m/s\ :sup:`2`. Defaults to Earth gravity.
    radius : float (optional)
        Planet radius in meters. Defaults to Earth's radius. 