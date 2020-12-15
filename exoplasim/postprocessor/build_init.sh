#!/bin/bash

cd res/src
tar -xvzf netcdf-cxx-legacy_4.2.orig.tar.gz
mv netcdf-cxx-4.2 legacy-netcdf-cxx-4.2
tar -xvzf netcdf-cxx4-4.2.tar.gz
cd ../
resdir=$(pwd)
cd src/netcdf-cxx4-4.2
export CPPFLAGS=$(nc-config --cflags)
export LDFLAGS=$(nc-config --libs)
export LD_LIBRARY_PATH=${LDFLAGS[@]:0:2}
./configure --prefix=$resdir
make
make install
cd ../legacy-netcdf-cxx-4.2
export CPPFLAGS=$(nc-config --cflags)
export LDFLAGS=$(nc-config --libs)
export LD_LIBRARY_PATH=${LDFLAGS[@]:0:2}
./configure --prefix=$resdir
make
make install
make check
cd $(dirname $resdir)
c++ -O2 -o burn7.x burn7.cpp -lm -Wl,-rpath=$resdir/lib -Lres/lib -lnetcdf_c++ -Ires/include