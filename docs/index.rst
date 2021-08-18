.. ExoPlaSim documentation master file, created by
   sphinx-quickstart on Mon Dec  7 19:50:40 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
ExoPlaSim Python API Documentation
==================================   

.. figure:: ../mixplanets.png
  :width: 100%
  :alt: Two rows of planets, progressing from yellow to blue from top left to bottom right. The top row appears to represent tidally-locked planets, while the bottom row appears to represent Earth-like planets.
  
  A range of planets modeled by ExoPlaSim, and postprocessed with SBDART. The top row consists of tidally-locked aquaplanets at T21 orbiting stars ranging from 2500 K to 4000 K, with orbital periods increasing with stellar mass. The bottom row consists of aquaplanets with 24-hour rotation at T42, orbiting stars ranging from 4000 K to 8000 K. 

.. _tut: tutorial.html

__ tut_

* :ref:`genindex`
* `Tutorial`__
* `Postprocessor`__
* :ref:`search`

.. _post: postprocessor.html
__ post_

Contents
========

.. toctree::
   :maxdepth: 3
   
   tutorial.rst
   postprocessor.rst
   source/exoplasim
..    :caption: Contents:


Created by Adiv Paradise

Copyright 2020, Distributed under the `General Public License`__. 

.. _GPL: LICENSE.html

__ GPL_

This API was written with Python 3 in mind, but should work with
Python 2 and outdated versions of NumPy. 

Requirements
------------

* Python (including development libraries, e.g. python-dev or python3.9-dev on Ubuntu--if using anaconda, these should already be included in your installation)
* numpy
* scipy (only needed for additional utilities)
* matplotlib (only needed for additional utilities)
* a Fortran compiler
* a C compiler
* (optionally) MPI libraries for those compilers

Optional Requirements
---------------------

* netCDF4 (for netCDF support)
* h5py (for HDF5 support)
* NetCDF-C Library (for legacy ``burn7`` support)
    
Installation
------------

::

    pip install exoplasim
    
OR::

    python setup.py install
    
    
If you know you will want to use NetCDF or HDF5 output formats,
you can install their dependencies at install-time:

::

    pip install exoplasim[HDF5]
    
OR::

    pip install exoplasim[netCDF4]
    
OR::

    pip install exoplasim[netCDF4,HDF5]
    
The first time you import the module and try to create a model
after either installing or updating, ExoPlaSim will run a 
configuration script, write the install directory into its 
source code, and compile the `pyfft` library.

.. burn7 NetCDF postprocessor. You must 
.. have NetCDF libraries available in the path when this happens.
.. The burn7 compilation process will build and compile a patched
.. version of the NetCDF libraries necessary for burn7--burn7 makes
.. use of features anachronistic to a particular version of NetCDF
.. that no longer exists.

You may also configure and compile the model manually if you wish
to not use the Python API, by entering the exoplasim/ directory
and running first configure.sh, then compile.sh (compilation flags
are shown by running ``./compile.sh -h``). 

.. The postprocessor and its
.. libraries can be compiled by entering ``exoplasim/postprocessor/`` and
.. running ``./build_init.sh``.

Most Common Error Modes
-----------------------

There are 3 major ways in which ExoPlaSim can crash. One is related
to installation, one is related to model compilation/configuration,
and one is related to numerical stability. 

If the most recent MOST_DIAG file appears to have run to completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. This can be one of **two** things: either the
.. netcdf libraries are not properly installed (the first time you tried
.. to create an ExoPlaSim model, this would have resulted in error output
.. when the burn7 postprocessor and its custom-patched netcdf library were
.. being built (the custom library depends on netcdf.h from a standard
.. netcdf library)), or you have output codes in the postprocessor namelist
.. that aren't supported. Check the codes first. If that's not it, you
.. can test whether or not it's a netcdf issue, by going to the run folder:
.. 
.. ::
.. 
..     ./burn7.x -n<example.nl>burnout MOST.00000 MOST.00000.nc
.. 
.. If the problem was with the netcdf libraries, you will get an error
.. complaining about a missing linked object. If this is the case,
.. install netcdf (or load the netcdf module if you're in a cluster environment),
.. and then rebuild the postprocessor as follows:::
.. 
..     burndir=$(python -c "import exoplasim; print(exoplasim.__path__)" | tail -c +3 | head -c -3)/postprocessor
..     cd $burndir
..     ./build_init.sh
.. 
.. This should (if it runs without errors) produce a new, functional ``burn7.x``.
.. You can test it by copying it to your run's working directory, and trying
.. the ``./burn7.x`` command given above once again.


If in the run folder, diagnostic files are produced that appear to have
made it all the way to the end of the year (there is a summary tag
giving time elapsed and that sort of thing), then the problem is likely
with the postprocessor. It is likely that the error output will be informative;
if it is not clear how to resolve, please let me (the developer) know. 

If the postprocessor itself is *not* the problem, then it's likely you somehow passed
incorrect output codes to the postprocessor. This is the most common scenario for
postprocessor-related crashes. Check your inputs for any errors. In particular, note 
that climatology outputs are not available if storm climatology was not enabled.

If ExoPlaSim crashed almost immediately without producing output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. py:module:: exoplasim.makestellarspec

If things crashed and burned immediately, it's likely a configuration
problem. Check to make sure you aren't using a restart file from a run
that used a different resolution, or stellar spectrum files that aren't
formatted correctly (use the :py:mod:`makestellarspec <exoplasim.makestellarspec>`
utility to format Phoenix spectra for ExoPlaSim), or boundary condition
``.sra`` files that aren't properly-formatted.

If ExoPlaSim ran for a while and then crashed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If things were fine until they weren't, then it's likely ExoPlaSim encountered
a numerical instability of some kind. Some of these are physical (e.g. you ran
a model at a thousand times Earth's insolation, and the oceans boiled, or the model
was too cold and the physics broke), while some are not (something happened to
violate the CFL condition for the given timestep, or an unphysical oscillation
wasn't damped properly by the dynamical core and it grew exponentially). If this
happens, either try a model configuration that is more physically reasonable,
or if the problem appears not to have been physical, try reducing the timestep
or increasing hyperdiffusion. Sometimes it also works to slightly adjust a model
parameter such as surface pressure by a fraction of a percent or less--just enough
to nudge the model out of whatever chaotic local minimum it ran into, but not
enough to qualitatively change the resulting climate. 

New in ExoPlaSim 3.0.0, there is a "crash-tolerant" run mode. With this mode enabled,
a runtime crash will result in rewinding 10 years and resuming. This deals with many
of the most frustrating problems related to numerical instability. However, due to
the potential for infinite loops, this is only recommended for advanced users.

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


A Note on NetCDF and the (deprecated) Burn7 Postprocessor
---------------------------------------------------------

As of ExoPlaSim 3.0.0, ``burn7`` is deprecated. It is included for
legacy support, but it has been replaced by ``pyburn`` for recommended
usage. ``pyburn`` has been designed as a drop-in replacement, so no
changes to existing codes and scripts are required. 

The Burn7 postprocessor requires the ``netcdfcpp.h`` header file. 
netcdf-cxx distributions later than version 4.2 no longer include
this file. A patched version of netcdf-cxx4-4.2 is shipped with
exoplasim, and will be built by default. However, doing so requires
that ``netcdf.h`` be available, typically from a netcdf C library.
If you encounter problems on first use, it is likely because the
C++ compilers can't find ``netcdf.h``, and you may need to adjust
the system path to include it. In a cluster environment, that may
involve a command such as ``module load netcdf``. On a personal Linux
computer, as long as netcdf is installed system-wide, this should
not be a problem. We have noted some issues building the patched 
netcdf library on newer versions of Mac OS X. The build process has
not been fully-tested on other verions of Mac or on Windows.


