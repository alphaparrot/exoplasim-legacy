==================
ExoPlaSim Tutorial
==================

In this tutorial, we will model the habitable zone terrestrial planet TOI 700 d, and take a look at some of the data. This tutorial assumes that you have installed ExoPlaSim successfully, and have matplotlib installed.

Setting Up 
==========

First thing's first: we want to import ExoPlaSim, and instantiate our :py:mod:`Model <exoplasim.Model>` instance. We want to create our model run in the folder "toi700d_run", and run it at T21 resolution on 4 CPUs.

>>> import exoplasim as exo
>>> toi700d = exo.Model(workdir="toi700d_run",modelname="TOI-700d",ncpus=4,resolution="T21")

If the appropriate executable does not yet exist, it will be compiled now. If this is the first time an ExoPlaSim model has been created, then a configuration script will be run first to locate the necessary compilers, and the NetCDF postprocessor will be built and compiled as well. We assign the model a descriptive name through the ``modelname`` argument, which is not strictly necessary, but will prove useful later.

Configuring the model for TOI-700d
----------------------------------

TOI 700 d was discovered by the TESS telescope in January 2020 `(Gilbert, et al 2020) <https://ui.adsabs.harvard.edu/link_gateway/2020AJ....160..116G/doi:10.3847/1538-3881/aba4b2>`_. It orbits TOI 700, a 3480 K M2V dwarf just over 100 lightyears away. TOI 700 has a luminosity of :math:`0.0233\pm0.0011` L\ :math:`_\odot`\ , and is relatively quiet. TOI 700 d has the following parameters:

+-----------------+-----------------------------------------------------------------------------------+
| Radius          | :math:`1.19\pm0.11` R\ :math:`_\oplus`                                            |
+-----------------+-----------------------------------------------------------------------------------+
| Mass            | :math:`1.72^{+1.29}_{-0.63}` M\ :math:`_\oplus`                                   |
+-----------------+-----------------------------------------------------------------------------------+
| Period          | :math:`37.4260^{+0.0007}_{-0.0010}` days                                          |
+-----------------+-----------------------------------------------------------------------------------+
| Semi-major Axis | :math:`0.163^{+0.0026}_{-0.0027}` AU                                              |
+-----------------+-----------------------------------------------------------------------------------+
| Incident Flux   | :math:`1367` W/m\ :math:`^2\left(\frac{L}{a^2}\right)\approx1199` W/m\ :math:`^2` |
+-----------------+-----------------------------------------------------------------------------------+
| Surface Gravity | :math:`9.81` m/s\ :math:`^2\left(\frac{M}{R^2}\right)\approx11.9` m/s\ :math:`^2` |
+-----------------+-----------------------------------------------------------------------------------+

We don't know anything else about the planet, so we'll have to make some assumptions about the atmosphere and surface. For simplicity, we'll assume that the surface is entirely ocean-covered, and that the atmospheric mass scales with planetary mass. We'll also assume that the atmosphere is N\ :sub:`2`\ , CO\ :sub:`2`\ , and H\ :sub:`2` \O. The surface pressure relative to Earth can therefore be estimated as follows:

.. math:: 
    p_s \approx \frac{g}{g_\oplus}\left(\frac{M}{M_\oplus}\right)\left(\frac{R_\oplus}{R}\right)^2
    
This gives a surface pressure of approximately 1.47 bars. With that figured out, we can proceed to configure the model (right now it is configured with the barest of defaults--you should always configure the model, even if you pass no non-default arguments).

>>> toi700d.configure(startemp=3480.0, flux=1167.0,                           # Stellar parameters
>>>                   eccentricity=0.,obliquity=0.,fixedorbit=True,           # Orbital parameters
>>>                   synchronous=True,rotationperiod=37.426,                 # Rotation
>>>                   radius=1.19,gravity=11.9,aquaplanet=True,               # Bulk properties
>>>                   pN2=1.47*(1-360e-6),pCO2=1.47*360e-6,ozone=False,       # Atmosphere
>>>                   timestep=30.0,snapshots=720,physicsfilter="gp|exp|sp")  # Model dynamics
>>> toi700d.exportcfg()

This command edits all the namelists and boundary condition files appropriately. The :py:func:`exportcfg() <exoplasim.Model.exportcfg>`command writes a portable text configuration file, by default named ``TOI-700d.cfg`` using the model's ``modelname`` parameter, that another user could use to replicate our model by simply running ``toi700d.loadconfig("TOI-700d.cfg")``.For a full description of the parameters we could have passed, see :py:func:`exoplasim.Model.configure() <exoplasim.Model.configure>`. Here is a brief overview of what each parameter did:

        startemp = 3480.0
            Specified the effective blackbody temperature of the star--in this case, 3480 K.
        flux = 1167.0
            Specified the incident flux (insolation or instellation) at the planet: 1167 W/m\ :math:`^2`
        eccentricity = 0.0
            We set the orbital eccentricity to 0.
        obliquity = 0.0
            We set the planet's axial tilt to 0.
        fixedorbit = True
            Here, we don't want the orbit precessing or anything, so we keep our orbit fixed.
        synchronous = True,
            By setting this flag, we have told ExoPlaSim that this is a tidally-locked model. The default is for the Sun to be fixed in place over 180° longitude.
        rotationperiod = 37.426
            Since the planet is tidally-locked, we assume its rotation period matches its orbital period, 37.426 days.
        radius = 1.19
            We set the planet's radius to 1.19 Earth radii.
        gravity = 11.9
            We set the surface gravity to 11.9 m/s\ :math:`^2`\ . Note that we do not specify the planet's mass directly, only the radius and surface gravity.
        aquaplanet = True
            Setting this flag deletes all surface boundary condition files and tells ExoPlaSim to initialize an ocean everywhere. The default is to have a mixed-layer depth of 50 meters.
        pN2 = 1.47*(1-360e-6)
            We want 1.47 bars **total**, but we want to include CO\ :sub:`2` as well. The surface pressure is the sum of the partial pressures, so we reduce pN\ :sub:`2` by the amount of CO\ :sub:`2` we want, the TOI 700 d equivalent of 360 :math:`\mu`\ bars. We could also skip the 1.47 scaling and set the pressure directly through its own argument.
        pCO2 = 1.47*360e-6
            We set the CO\ :sub:`2` partial pressure to its Earth level in bars, scaled up.
        ozone = False
            Since we are not assuming an oxygenated atmosphere (and some studies dispute how much ozone could be produced from an oxygenated atmosphere around an M dwarf anyway), we assume there will be no forcing from ozone. Tidally-locked models in ExoPlaSim are more stable without ozone anyway.
        timestep = 30.0
            Tidally-locked climates are stlightly more extreme than Earth-like climates, so rather than the default 45-minute timestep, we use 30 minutes.
        snapshots = 720
            Here we tell ExoPlaSim to write snapshot outputs every 720 timesteps (15 days). These snapshots show us the climate at a particular instant in time, and are therefore necessary for any observational postprocessing (any time-integrated observation is an average of photons that passed through the atmosphere as it was for a brief moment, not through the time-averaged atmosphere--this is mainly important for clouds). It's usually a good idea to write a snapshot every 15 days (twice a month), so scale based on the timestep. The default is to write every 480 timesteps, which is 15 days when a timestep is 15 minutes.
        physicsfilter = "gp|exp|sp"
            Tidally-locked models can be subject to large-scale Gibbs oscillations on the night side, due to the strong dipole moment of the forcing and axial symmetry of the iceline. **All models will struggle to reproduce sharp features accurately**. ExoPlaSim merely struggles in an extremely visible way. Fortunately, we can mitigate this to an acceptable level with the use of *physics filters*. These are mathematical filters included in the dynamical core at the spectral transform stage. Here we have told ExoPlaSim to use an exponential filter, and to apply it both at the transform from gridpoint space to spectral space, and at the transform from spectral space back to gridpoint space. For more details on the choice of filter and how they work, see :py:func:`exoplasim.Model.configure() <exoplasim.Model.configure>`. For Earth-like models that aren't tidally-locked, physics filters are usually not necessary.
            
Running the Model
-----------------

Now that we have configured the model, it's time to run it! This demo is intended to be something you can run on your laptop (thus specifying only 4 CPUs), so to make sure you have something to look at when you come back from your lunch break, let's just run for 10 years. On my laptop with 4 cores, a year takes just over 6 minutes. Note that on HPC architecture with 16 cores, a year often takes less than a minute.

>>> toi700d.run(years=10,crashifbroken=True)

The ``crashifbroken`` flag simply means that if something goes wrong, the model will crash in a slightly cleaner, Pythonic way. Note that a problem with the postprocessor will get flagged as a crash just like an actual model crash--in most cases, the model is salvageable if you figure out what went wrong with the postprocessor.
            
    