.. -*- coding:utf-8 -*-

==================================
ExoPlaSim-Legacy Python API README
==================================

Created by Adiv Paradise

Copyright 2020, Distributed under the General Public License

This API was written with Python 3 in mind, but should work with
Python 2 and outdated versions of NumPy. 

Read the full documentation at http://exoplasim.readthedocs.io.

Requirements
------------
    
*   numpy
*   scipy
*   matplotlib (only needed for additional utilities)
*   GNU C (gcc/g++) and Fortran (gfortran) compilers (for Python utilities)
*   (optionally) Other compilers whose use you prefer for the model itself
*   (optionally) MPI libraries for those compilers
*   netCDF4 (optional)
*   h5py (optional)
    
Installation
------------

::

    pip install exoplasim-legacy
    
OR::

    python setup.py install
    
The first time you import the module and try to create a model
after either installing or updating, ExoPlaSim will run a 
configuration script, write the install directory into its 
source code, and (if applicable) compile the burn7 NetCDF postprocessor.

Multiple output formats are supported by the built-in `pyburn`
postprocessor. If you wish to use HDF5 or NetCDF output formats, you
will need the netCDF4-python and h5py libraries, respectively. You
can ensure these are included at install-time by specifying them:

::
    pip install exoplasim-legacy[netCDF4]
    
OR::
    pip install exoplasim-legacy[HDF5]
    
OR::
    pip install exoplasim-legacy[netCDF4,HDF5]

You may also configure and compile the model manually if you wish
to not use the Python API, by entering the exoplasim/ directory
and running first configure.sh, then compile.sh (compilation flags
are shown by running ``./compile.sh -h``). The postprocessor and its
libraries can be compiled by entering ``exoplasimlegacy/postprocessor/`` and
running ``./build_init.sh``.

burn7 compilation
-----------------
You must have NetCDF libraries available in the path to build burn7.
The burn7 compilation process will build and compile a patched
version of the NetCDF libraries necessary for burn7--burn7 makes
use of features anachronistic to a particular version of NetCDF
that no longer exists.


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

    import exoplasimlegacy as exo
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
