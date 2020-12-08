.. -*- coding:utf-8 -*-

===========================
ExoPlaSim Python API README
===========================

Created by Adiv Paradise

Copyright 2020, Distributed under the General Public License

This API was written with Python 3 in mind, but should work with
Python 2 and outdated versions of NumPy. 

Requirements
------------
    
* netCDF4
* numpy
* a Fortran compiler
* a C compiler
* (optionally) MPI libraries for those compilers
    
Installation
------------

::

    pip install exoplasim
    
OR::

    python setup.py install
    
The first time you import the module and try to create a model
after either installing or updating, ExoPlaSim will run a 
configuration script, write the install directory into its 
source code, and compile the burn7 NetCDF postprocessor. You must 
have NetCDF libraries available in the path when this happens.
The burn7 compilation process will build and compile a patched
version of the NetCDF libraries necessary for burn7--burn7 makes
use of features anachronistic to a particular version of NetCDF
that no longer exists.

You may also configure and compile the model manually if you wish
to not use the Python API, by entering the exoplasim/ directory
and running first configure.sh, then compile.sh (compilation flags
are shown by running ``./compile.sh -h``). The postprocessor and its
libraries can be compiled by entering ``exoplasim/postprocessor/`` and
running ``./build_init.sh``.

PlaSim Documentation
--------------------

Original PlaSim documentation is available in the exoplasim/docs/
folder.

Usage
-----

To use the ExoPlaSim Python API, you must import the module, create
a Model or one of its subclasses, call its configure method and/or
modify method, and then run it. 

Basic example:::

    import exoplasim as exo
    mymodel = exo.Model(workdir="mymodel_testrun",modelname="mymodel",resolution="T21",layers=10,ncpus=8)
    mymodel.configure()
    mymodel.exportcfg()
    mymodel.run(years=100,crashifbroken=True)
    mymodel.finalize("mymodel_output")
    
In this example, we initialize a model that will run in the directory
"mymodel_testrun", and has the name "mymodel", which will be used to
label output and error logs. The model has T21 resolution, or 32x64,
10 layers, and will run on 8 CPUs. By default, the compiler will use
8-byte precision. 4-byte may run slightly faster, but possibly at the
cost of reduced stability. If there are machine-specific optimization
flags you would like to use when compiling, you may specify them as a
string to the optimization argument, e.g. ``optimization='mavx'``. ExoPlaSim
will check to see if an appropriate executable has already been created,
and if not (or if flags indicating special compiler behavior such as 
debug=True or an optimization flag are set) it will compile one. We then
configure the model with all the default parameter choices, which means
we will get a model of Earth. We then export the model configurations
to a ``.cfg`` file (named automatically after the model), which will allow
the model configuration to be recreated exactly by other users. We 
run the model for 100 years, with error-handling enabled. Finally, we 
tell the model to clean up after itself. It will take the most recent 
output files and rename them after the model name we chose, and delete 
all the intermediate output and configuration files. 

Objects
-------

``Model(resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,recompile=False,optimization=None,mars=False,workdir="most",source=None,modelname="MOST_EXP")``
    Initialize an ExoPlaSim model in a particular directory. If the necessary executable does not yet exist, compile it.
      
:resolution: 
    The resolution of the model. Options are T21, T42, T63, T85, 
    T106, T127, and T170, corresponding to 32, 64, 96, 128, 160, 
    192, and 256 latitudes respectively, and twice as many 
    longitudes. ExoPlaSim has been tested and validated most 
    extensively at T21 and T42. Higher resolutions will take 
    considerable time to run.
                  
:layers: 
    The number of vertical layers in the model atmosphere. The default
    is 10, but PlaSim has been used with 5 layers in many studies.
    More layers are supported, but not recommended except at higher
    resolutions.
              
:ncpus: 
    The number of MPI processes to use, typically the number of cores
    available. If ncpus=1, MPI will not be used.
       
:precision: 
    Either 4 or 8--specifies the number of bytes for a Fortran real.

:debug: 
    Either True or False. If True, compiler optimizations are disabled
    and the code is compiled with debugging flags enabled that will
    allow line-by-line tracebacks if ExoPlaSim crashes. Only use for
    development purposes.
       
:inityear: 
    The number to use for the initial model year (default 0).

:recompile: 
    True/False flag used to force a recompile. Cannot force the 
    model to skip compilation if the executable does not exist or
    compilation-inducing flags are set.
           
:optimization: 
    Fortran compiler arguments for optimization. ANY compiler
    flags can be passed here, but it's intended for optimization
    flags. Setting this will trigger a recompile.
              
:mars: 
    True/False. If True, will use Mars-specific routines.

:workdir: 
    The directory in which to construct the model.

:source: 
    The directory in which to look for executables, namelists, 
    boundary conditions, etc. If not set, will default to exoplasim/plasim/run/.
        
:modelname: 
    The name to use for the model and its output files when finished.
      
``TLmodel(resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,recompile=False,optimization=None,mars=False,workdir="most",source=None,modelname="MOST_EXP")``
    Identical to Model, except configuration options suitable for
    tidally-locked models are the default when configure() is called.
                
``TLaquaplanet(resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,recompile=False,optimization=None,mars=False,workdir="most",source=None,modelname="MOST_EXP")``
    Identical to Model, except configuration options suitable for
    tidally-locked models are the default when configure() is called,
    and the surface is entirely ocean-covered. Specifically, a 30-minute
    timestep, snapshot outputs every 720 timesteps, eccentricity=0.0,
    0-degree obliquity, exponential physics filtering, fixed orbital
    parameters, and no ozone. All these defaults can be overridden.
    
``TLlandplanet(resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,recompile=False,optimization=None,mars=False,workdir="most",source=None,modelname="MOST_EXP")``
    Identical to Model, except configuration options suitable for
    tidally-locked models are the default when configure() is called,
    and the surface is entirely land-covered. Specifically, a 30-minute
    timestep, snapshot outputs every 720 timesteps, eccentricity=0.0,
    0-degree obliquity, exponential physics filtering, fixed orbital
    parameters, and no ozone. All these defaults can be overridden.
    
``Earthlike(resolution="T21",layers=10,ncpus=4,precision=8,debug=False,inityear=0,recompile=False,optimization=None,mars=False,workdir="most",source=None,modelname="MOST_EXP")``
    Identical to Model, except configuration options common for
    Earth-like models are the default when configure is called--
    specifically, 45-minute timestep, snapshot output reporting 
    every 480 timesteps, and a model top pinned to 50 mbar. All 
    these defaults can be overridden.

Methods
-------

``.configure(noutput=True,flux=1367.0,startemp=None,starspec=None,pH2=None,pHe=None,pN2=None,pO2=None,pCO2=None,pAr=None,pNe=None,pKr=None,pH2O=None,gascon=None,pressure=None,pressurebroaden=True,vtype=0,rotationperiod=24.0,synchronous=False,substellarlon=180.0,year=None,glaciers={"toggle":False,"mindepth":2.0,"initialh":-1.0},restartfile=None,gravity=9.80665,radius=1.0,eccentricity=None,obliquity=None,lonvernaleq=None,fixedorbit=False,orography=None,seaice=True,co2weathering=False,evolveco2=False,physicsfilter=None,filterkappa=8.0,filterpower=8,filterLHN0=15.0,diffusionwaven=None,qdiffusion=None,tdiffusion=None,zdiffusion=None,ddiffusion=None,diffusionpower=None,erosionsupplylimit=None,outgassing=50.0,snowicealbedo=None,twobandalbedo=False,maxsnow=None,soilalbedo=None,oceanalbedo=None,oceanzenith="ECHAM-3",wetsoil=False,soilwatercap=None,aquaplanet=False,desertplanet=False,soilsaturation=None,drycore=False,ozone=True,cpsoil=None,soildepth=1.0,mldepth=50.0,writefrequency=None,modeltop=None,stratosphere=False,tropopause=None,timestep=45.0,runscript=None,columnmode=None,highcadence={"toggle":0,"start":320,"end":576,"interval":4},snapshots=None,resources=[],landmap=None,stormclim=False,nstorms=4,stormcapture={"VITHRESH":0.145,"GPITHRESH":0.37,"VMXTHRESH":33.0,"LAVTHRESH":1.2e-5,"VRMTHRESH":0.577,"MINSURFTEMP":298.15,"MAXSURFTEMP":373.15,"WINDTHRESH":33.0,"SWINDTHRESH":20.5,"SIZETHRESH":30,"ENDTHRESH":16,"MINSTORMLEN":256,"MAXSTORMLEN":1024,"NKTRIGGER":0,"toggle":0},topomap=None,threshold=5.0e-4,otherargs={})``
    Configure model boundary conditions and namelist files. Defaults are appropriate
    for an Earth model. The various options will be organized below by category.
    
Model Operation
~~~~~~~~~~~~~~~

    :noutput: 
        True/False. Whether or not model output should be written.
        restartfile: Path to a restart file to use for initial conditions. Can be None.
    
    :writefrequency: 
        How many times per day ExoPlaSim should write output. Ignored by
        default--default is to write time-averaged output once every 5 days.
        
    :timestep: 
        Model timestep. Defaults to 45 minutes.
        
    :runscript: 
        A Python function that accepts a Model object as its first argument. This
        is the routine that will be run when you issue the Model.run() command.
        Any keyword arguments passed to run() will be forwarded to the specified
        function. If not set, the default internal routine will be used.
        
    :snapshots: 
        How many timesteps should elapse between snapshot outputs. If not set,
        no snapshots will be written.
        
    :highcadence: 
        A dictionary containing the following arguments:
         ``toggle``:    1/0. Whether or not high-cadence output should be written (1=yes).
         ``start``:     Timestep at which high-cadence output should begin.
         ``end``:       Timestep at which high-cadence output should end.
         ``interval``:  How many timesteps should elapse between high-cadence outputs.
        
    :threshold: 
        Energy balance threshold model should run to, if using runtobalance().
        Default is <0.05 W/m^2/yr average drift in TOA and surface energy balance
        over 45-year timescales.
            
    :resources: 
        A list of paths to any additional files that should be available in the
        run directory.
            
    :otherargs: 
        Any namelist parameters not included by default in the configuration options.
        These should be passed as a dictionary, with "PARAMETER@namelist" as the
        form of the dictionary key, and the parameter value passed as a string.
        e.g. ``otherargs={"N_RUN_MONTHS@plasim_namelist":'4',"NGUI@plasim_namelist:'1'}``
           
Model Dynamics
~~~~~~~~~~~~~~

    :columnmode: 
        Can be "-", "clear", "static", "static|clear", or "clear|static". The 
                inclusion of 'static' will disable horizontal advection, forcing ExoPlaSim
                into a column-only mode of operation. The inclusion of 'clear' will disable
                the radiative effects of clouds.
                
    :drycore: 
        True/False. If True, evaporation is turned off, and a dry atmosphere will
        be used.
             
    :physicsfilter: 
        If not an empty string, specifies the physics filter(s) to be used. Filters
        can be used during the transform from gridpoint to spectral (``"gp"``), and/or
        during the transform from spectral to gridpoint (``"sp"``). Filter types are:
        
             ``"none"``:   No filter
             ``"cesaro"``: f(n) = 1 - n/(N+1)           (Cesaro filter)
             ``"exp"``:    f(n) = exp(-k*(n/N)^y)          (Exponential filter)
             ``"lh"``:     f(n) = exp(-[n(n+1)/(g(g+1))]^2)   (Lander-Hoskins filter)
             
        Where ``n`` is the wavenumber, ``N`` is the truncation wavenumber (e.g. 21 for
        T21), ``k`` is ``filterkappa``, ``y`` is ``filterpower``, and ``g`` is ``filterLHN0``.
        Combinations of filter types and times should be combined with a ``|``,
        e.g. ``physicsfilter="gp|exp|sp"`` or ``physicsfilter="gp|cesaro"``.
                    
    :filterkappa: 
        A constant to be used with the exponential filter. Default is 8.0.
    
    :filterpower: 
        A constant integer to be used with the exponential filter. Default is 8.
    
    :filterLHN0: 
        The constant used in the denominator of the Lander-Hoskins Filter. Default
        is 15; typically chosen so f(N)=0.1.
                
    :diffusionwaven: 
        The critical wavenumber beyond which hyperdiffusion is applied. Default
        is 15 for T21.
                
    :qdiffusion: 
        Timescale for humidity hyperdiffusion in days. Default for T21 is 0.1.
    
    :tdiffusion: 
        Timescale for temperature hyperdiffusion in days. Default for T21 is 5.6.
    
    :zdiffusion: 
        Timescale for vorticity hyperdiffusion in days. Default for T21 is 1.1.
    
    :ddiffusion: 
        Timescale for divergence hyperdiffusion in days.. Default for T21 is 0.2.
    
    :diffusionpower: 
        integer exponent used in hyperdiffusion. Default is 2 for T21.
        
Radiation
~~~~~~~~~

    :flux: 
        Incident stellar flux in W/m^2. Default 1367 for Earth.
    
    :startemp: 
        Effective blackbody temperature for the star. Not used if not set.
    
    :starspec: 
        Spectral file for the stellar spectrum. Should have two columns and 965 rows,
        with wavelength in the first column and radiance or intensity in the second.
        A similarly-named file with the "_hr.dat" suffix must also exist and have 
        2048 wavelengths.
              
    :twobandalbedo: 
        True/False. If True, separate albedos will be calculated for each of the
        two shortwave bands. If False (default), a single broadband albedo will be
        computed and used for both.
                   
    :synchronous: 
        True/False. If True, the Sun is fixed to one longitude in the sky.
    
    :substellarlon: 
        The longitude of the substellar point, if synchronous==True. Default 180°
    
    :pressurebroaden: 
        True/False. If False, pressure-broadening of absorbers no longer depends
        on surface pressure. Default is True
                     
    :ozone: 
        True/False. Whether or not forcing from stratospheric ozone should be included.
    
    :snowicealbedo: 
        A uniform albedo to use for all snow and ice.
    
    :soilalbedo: 
        A uniform albedo to use for all land.
    
    :wetsoil: 
        True/False. If True, land albedo depends on soil moisture (wet=darker).
    
    :oceanalbedo: 
        A uniform albedo to use for the ocean.
    
    :oceanzenith: 
        The zenith-angle dependence to use for blue-light reflectance from the ocean.
        Can be ``'Lambertian'``/``'uniform'``, ``'ECHAM-3'``/``'plasim'``/``'default'``, or ``'ECHAM-6'``.
        The default is ``'ECHAM-3'`` (synonymous with ``'plasim'`` and ``'default'``), which is
        the dependence used in the ECHAM-3 model.
                     
Orbital Parameters
~~~~~~~~~~~~~~~~~~

     :year: 
        Number of 24-hour days in a sidereal year. Not necessary if eccentricity and 
        obliquity are zero. Defaults if not set to ~365.25 days
           
     :rotationperiod: 
        Planetary rotation period, in days. Default is 1.0.
     
     :eccentricity: 
        Orbital eccentricity. If not set, defaults to Earth's (0.016715)
     
     :obliquity: 
        Axial tilt, in degrees. If not set, defaults to Earth's obliquity (23.441°).
     
     :lonvernaleq: 
        Longitude of periapse, measured from vernal equinox, in degrees. If 
        not set, defaults to Earth's (102.7°).
                  
     :fixedorbit: 
        True/False. If True, orbital parameters do not vary over time. If False,
        variations such as Milankovich cycles will be computed by PlaSim.
        
Planet Parameters
~~~~~~~~~~~~~~~~~

     :gravity: 
        Surface gravity, in m/s^2. Defaults to 9.80665 m/s^2.
     
     :radius: 
        Planet radius in Earth radii. Default is 1.0.
     
     :orography: 
        If set, a scaling factor for topographic relief. If ``orography=0.0``, topography
        will be zeroed-out.
               
     :aquaplanet: 
        True/False. If True, the surface will be entirely ocean-covered.
     
     :desertplanet: 
        True/False. If True, the surface will be entirely land-covered.
     
     :seaice: 
        True/False. If False, disables radiative effects of sea ice (although sea ice 
        itself is still computed).
             
     :landmap: 
        Path to a ``.sra`` file containing a land mask for the chosen resolution.
     
     :topomap: 
        Path to a ``.sra`` file containing geopotential height map. Must include landmap.
        
Atmosphere
~~~~~~~~~~

     :gascon: 
         Effective gas constant. Defaults to 287.0 (Earth), or the gas constant
         corresponding to the composition specified by partial pressures.
             
     :vtype: 
         Type of vertical discretization. Can be:
         0   Pseudolinear scaling with pressure that maintains resolution near the ground.
         1   Linear scaling with pressure.
         2   Logarithmic scaling with pressure (resolves high altitudes)
         3   Pseudologarithmic scaling with pressure that preserves resolution near the ground.
         4   Pseudolinear scaling with pressure, pinned to a specified top pressure.
         5   If >10 layers, bottom 10 as if ``vtype=4``, and upper layers as if ``vtype=2``.
         
     :modeltop: 
         Pressure of the top layer
     
     :tropopause: 
         If stratosphere is being included, pressure of the 10th layer (where scheme
         switches from linear to logarithmic).
                 
     :stratosphere: 
         True/False. If True, vtype=5 is used, and model is discretized to include
         a stratosphere.
         
     :pressure: 
            Surface pressure in bars, if not specified through partial pressures.
        
Gas Partial Pressures
^^^^^^^^^^^^^^^^^^^^^

    Partial pressures of individual gases can be specified. If pressure and gascon
    are not explicitly set, these will determine surface pressure, mean molecular
    weight, and effective gas constant. Note however that Rayleigh scattering assumes
    an Earth-like composition, and the only absorbers explicitly included in the 
    radiation scheme are CO2 and H2O.
    
     :pH2:   
        H2 partial pressure in bars.
     
     :pHe:   
        He partial pressure in bars.
     
     :pN2:  
        N2 partial pressure in bars.
     
     :pO2:  
        O2 partial pressure in bars.
     
     :pH2:  
        H2 partial pressure in bars.
     
     :pAr:  
        Ar partial pressure in bars.
     
     :pNe:  
        Ne partial pressure in bars.
     
     :pKr:  
        Kr partial pressure in bars.
     
     :pCO2:  
        CO2 partial pressure in bars. This gets translated into a ppmv concentration, so if you want to specify/vary CO2 but don't need the other gases, specifying pCO2, pressure, and gascon will do the trick. In most use cases, however, just specifying pN2 and pCO2 will give good enough behavior.
            
    :pH2O:  
        H2O partial pressure in bars. This is only useful in setting the gas constant and surface pressure; it will have no effect on actual moist processes.
        
Surface Parameters
~~~~~~~~~~~~~~~~~~

    :mldepth: 
        Depth of the mixed-layer ocean. Default is 50 meters.
    
    :soildepth: 
        Scaling factor for the depth of soil layers (default total of 12.4 meters)
    
    :cpsoil: 
        Heat capacity of the soil, in J/m^3/K. Default is 2.4*10^6.
    
    :soilwatercap: 
        Water capacity of the soil, in meters. Defaults to 0.5 meters
    
    :soilsaturation: 
        Initial fractional saturation of the soil. Default is 0.0 (dry).
    
    :maxsnow: 
        Maximum snow depth (Default is 5 meters; set to -1 to have no limit).
        
Additional Physics
~~~~~~~~~~~~~~~~~~

    :co2weathering: 
        True/False. Toggles whether or not carbon-silicate weathering should be
        computed. Default is False.
        
    :evolveco2: 
        True/False. If co2weathering==True, toggles whether or not the CO2 partial
        pressure should be updated every year. Usually the change in pCO2 will be 
        extremely small, so this is not necessary, and weathering experiments try
        to estimate the average weathering rate for a given climate in order to 
        interpolate timescales between climates, rather than modelling changes in CO2
        over time directly.
        
    :outgassing: 
        The assumed CO2 outgassing rate in units of Earth outgassing. Default is 1.0.
        
    :erosionsupplylimit: 
        If set, the maximum CO2 weathering rate per year permitted by
        erosion, in ubars/year. This is not simply a hard cutoff, but follows
        Foley 2015 so high weathering below the cutoff is also reduced.
    
    :glaciers: 
        A dictionary containing the following arguments:
        
        toggle:  True/False. Whether or not glaciers should be allowed to grow or shrink in thickness, or be formed from persistent snow on land.
        mindepth:  The minimum snow depth in meters of liquid water equivalent that must persist year-round before the grid cell is considered glaciated. Default is 2 meters.
        initialh:  If >=0, covers the land surface with ice sheets of a height given in meterss. If -1, no initial ice sheets are assumed.

                  
    :stormclim: 
        True/False. Toggles whether or not storm climatology (convective available
        potential energy, maximum potential intensity, ventilation index, etc)
        should be computed. If True, output fields related to storm climatology 
        will be added to standard output files. Enabling this mode currently roughly
        doubles the computational cost of the model. This may improve in future 
        updates. Refer to Paradise, et al 2021 for implementation description. 
        
    :stormcapture: 
        A dictionary containing arguments controlling when high-cadence output
        is triggered by storm activity. This dictionary must contain 'toggle', which
        can be either 1 or 0 (yes or no). It may also contain any namelist
        parameters accepted by hurricanemod.f90, including the following:
        
        NKTRIGGER:  0/1 (no/yes). Whether or not to use the Komacek, et al 2020 conditions for hurricane cyclogenesis as the output trigger. Default is no.
        VITHRESH:   (nktrigger) Ventilation index threshold for nktrigger output. Default 0.145
        VMXTHRESH:  (nktrigger) Max potential intensity threshold for nktrigger output.Default 33 m/s
        LAVTHRESH:  (nktrigger) Lower-atmosphere vorticity threshold for nktrigger output. Default 1.2*10^-5 s^-1
        VRMTHRESH:  (unused) Ventilation-reduced maximum intensity threshold. Default 0.577
        GPITHRESH:  (default) Genesis Potential Index threshold. Default 0.37.
        MINSURFTEMP:  (default) Min. surface temperature for storm activity. Default 25C
        MAXSURFTEMP:  (default) Max. surface temperature for storm activity. Default 100C
        WINDTHRESH:   (default) Lower-atmosphere maximum wind threshold for storm activity.  Default 33 m/s
        SWINDTHRESH:  (default) Minimum surface windspeed for storm activity. Default 20.5 m/s
        SIZETHRESH:  (default) Minimum number of cells that must trigger to start outputDefault 30
        ENDTHRESH:  (default) Minimum number of cells at which point storm output ends.Default 16
        MINSTORMLEN:  (default) Minimum number of timesteps to write output. Default 256
        MAXSTORMLEN:  (default) Maximum number of timesteps to write output. Default 1024
    
    Note that actual number of writes will be stormlen/interval, as set in
    highcadence. This interval defaults to 4, so 64 writes minimum, 256 max

Other Model Methods
-------------------
 
``.modify(**kwargs)``   
    Modify any already-configured parameters. All parameters accepted by configure() can
    be passed as arguments.
    
``.exportcfg(filename=None)``
    Write the current model configuration to a text file. This file can be shared and used by
    other users to recreate your model configuration.
    
        :filename: 
            Path to the file that should be written. If None (default), <modelname>.cfg
            will be created in the working directory.
                  
``.loadconfig(configfile)``
    Load a previously-exported configuration file and configure the model accordingly.
    
        :configfile: 
            Path to the configuration file to load
        
``.save(filename=None)``
    Save the current Model object to a NumPy save file. The model object can then be 
    reinstantiated using ``numpy.load(savefile).item()``. Note that these are often not
    portable between versions of Python or machine architectures, so their use is only
    recommended internally. For sharing with other users, it is recommended that you use
    the ``.exportcfg()`` function.
    
``.run(\*\*kwargs)``  
   Run the model's run routine. This may have been passed as runscript when the model was
   created, or it could be the model's internal ._run() routine. That method is described
   below:
    
.. code:: python    

        ._run(years=1,postprocess=True,crashifbroken=False,clean=True)
        
            :years: 
                Number of years to run
                
            :postprocess: 
                True/False. Whether or not NetCDF files should be produced on-the-fly
                
            :crashifbroken: 
                True/False. If True, use Pythonic error handling
                
            :clean: 
                True/False. If True, delete raw output files once NetCDF files are made

``.runtobalance(threshold=None,baseline=50,maxyears=300,minyears=75,timelimit=None,crashifbroken=True,clean=True)``
    Run the model until energy balance equilibrium is reached at the top and surface.
        
        :threshold: 
            If specified, overrides the threshold set by ``.config()``. The model will run
            until the energy balance at the top and surface drifts by less than this
            amount per year over a given baseline.
                   
        :baseline: 
            The number of years over which to evaluate energy balance drift. Default 50
        
        :maxyears: 
            The maximum number of years to run before returning. Default 300. This is
            useful if you are running on a scratch disk with limited space.
                  
        :minyears: 
            The minimum number of years to run before determining that the model is in
            equilibrium.
                  
        :timelimit: 
            If set, maxyears will be revised each year based on the average minutes
            per year thus far, to try to avoid going over the time limit, which should
            be given in minutes.
                   
        :crashifbroken: 
            True/False. If True, Pythonic error handling is enabled. Default True.
        
        :clean: 
            True/False. If True, raw output is deleted once postprocessed. Default True.
        
``.postprocess(inputfile,namelist,log="postprocess.log",crashifbroken=False)``
    Produce NetCDF output from an input file, using a specified postprocessing namelist. 
    
        :inputfile: 
            The raw output file to be processed
            
        :namelist: 
            The burn7 namelist to use
            
        :log: 
            The log file to which burn7 should output standard output and errors
            
        :crashifbroken: 
            True/False. If True, exoplasim will run .integritycheck() on the file.
        
``.integritycheck(ncfile)``
    Check an output file to see it contains the expected variables and isn't full of NaNs.
    If the file does not exist, exoplasim will attempt to create it using the postprocessor.
    If the file does not have the expected variables or is full of trash, an exception will
    be raised. If the file is fine, this function returns a 1. If the file did not exist and
    cannot be created, this function will return a 0. 
    
        :ncfile: 
            The output file to check.
    
``.finalize(outputdir,allyears=False,keeprestarts=False,clean=True)``
    Move outputs and optionally restarts to a specified output directory. If more than the final
    year of output is being kept, a folder will be created in the output directory using the 
    model name. Otherwise, finalized files will be renamed using the model name.
    
        :outputdir: 
            Directory in which to put output.
            
        :allyears: 
            True/False. If True, output from all years will be kept, in a directory in
            outputdir named with the model name. Otherwise, the most recent year will be
            kept in outputdir, using the model name. Default False.
            
        :keeprestarts: 
            True/False: If True, restart files will be kept as well as output files.
            Default False.
            
        :clean: 
            True/False. If True, the original working directory will be deleted after files
            are moved. Default True.
        
``.gethistory(key="ts",mean=True,layer=-1)``
    Return an array of the global annual average for the provided output key, for each year
    of output.
    
        :key: 
            One of the supported output variables.
            
        :mean: 
            True/False. If False, compute the global sum instead of the global mean.
            
        :layer: 
            If the output field specified has 3 spatial dimensions, this is the level
            that will be returned. Level indexing follows Python rules, so negative
            numbers are allowed. Item 0 is the top of the atmosphere; -1 is the bottom.

``.get(year,snapshot=False,highcadence=False)``
    Return an open NetCDF data object for the given year. Defaults is to return time-averaged
    output.
    
        :year: 
            Integer year of output to return
            
        :snapshot: 
            True/False. If True, return the snapshot version.
            
        :highcadence: 
            True/False. If True, return the high-cadence version.
    
``.inspect(variable,year=-1,ignoreNaNs=True,snapshot=False,highcadence=False,savg=False,tavg=False,layer=None)``
    Return a given output variable from a given year, with optional averaging parameters.
    
        :variable: 
            The name of the variable to return.
            
        :year: 
            Which year of output to return. Year indexing follows Pythonic rules. If the model
            has been finalized, only the final year of output will be returned.
            
        :ignoreNaNs: 
            True/False. If True, use NaN-tolerant numpy functions.
            
        :snapshot: 
            True/False. If True, use snapshot output instead of time-averaged.
            
        :highcadence: 
            True/False. If True, use high-cadednce output instead of time-averaged.
            
        :savg: 
            True/False. If True, compute the spatial average. Default False
            
        :tavg: 
            True/False. If True, compute the annual average. Default False
            
        :layer: 
            If specified and data has 3 spatial dimensions, extract the specified layer. If
            unspecified and data has 3 spatial dimensions, the vertical dimension will be
            preserved (even if spatial averages are being computed).
