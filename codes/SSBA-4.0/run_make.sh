#!/bin/bash

clear

rm -rf build

mkdir build && cd build && cmake .. && make

clear

cd Apps

DATA_FILE=~/Dataset/problem-308-195089-pre.txt
PARAM_FILE=~/BundleAdjustment/codes/SSBA-4.0/params.txt


#USAGE:
# ./bundle_large_lifted_schur <sparse reconstruction file> <parameter file>
# params.txt : MAX_GMRES_ITERATIONS RESTARTS_FOR_GMRES TOLERANCE
./bundle_large_lifted_schur  $DATA_FILE $PARAM_FILE
#./bundle_large_lifted_schur ~/Dataset/problem-49-7776-pre.txt
#./bundle_large ~/Dataset/problem-49-7776-pre.txt
