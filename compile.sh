#!/bin/bash

prec=4
resolution="t21"
levels=10
ncpus=4
debug=""
optimization=""

while getopts "p:r:v:n:O:d" opt; do
    case $opt in
        p)
            case $OPTARG in
                4)
                    prec=4 ;;
                8)
                    prec=8 ;;
                single)
                    prec=4 ;;
                double)
                    prec=8 ;;
                \?)
                    echo "INVALID PRECISION. Reverting to single precision."
                    ;;
            esac
            ;;
        r)
            case $OPTARG in
                T21)
                    resolution="t21"
                    ;;
                T42)
                    resolution="t42"
                    ;;
                T63)
                    resolution="t63"
                    ;;
                32)
                    resolution="t21"
                    ;;
                64)
                    resolution="t42"
                    ;;
                96)
                    resolution="t63"
                    ;;
                \?)
                    resolution="t21"
                    echo "INVALID RESOLUTION PASSED! Reverting to T21."
                    ;;
            esac ;;
        v)
            levels=$OPTARG
            ;;
        n)
            ncpus=$OPTARG
            ;;
        d)
            debug="-ggdb -traceback -fpe0 "
            ;;
        O)
            optimization="-"$OPTARG
            ;;
        \?)
            echo "UNRECOGNIZED ARGUMENT PASSED: "$OPTARG
            ;;
    esac
done

echo "PRODUCING: "$optimization" -r"$prec" "$debug"-o most_plasim_"$resolution"_l"$levels"_p"$ncpus".x"
executable="most_plasim_"$resolution"_l"$levels"_p"$ncpus".x"

cd plasim/bld/
cp -p ../src/* .
[ ! -e MPI ] && rm -f *.o *.mod *.x
touch MPI
cat ../../most_compiler_mpi ../../most_debug_options make_plasim > makefile
make -e
./most_snow_build4
./most_ice_build4
cp plasim.x ../bin/$executable