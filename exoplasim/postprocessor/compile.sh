#!/bin/bash

resdir=$(pwd)/res
c++ -O2 -o burn7.x burn7.cpp -lm -Wl,-rpath=$resdir/lib -Lres/lib -lnetcdf_c++ -Ires/include