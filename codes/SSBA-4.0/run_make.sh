#!/bin/bash

clear

rm -rf build

mkdir build && cd build && cmake .. && make

clear

cd Apps


#USAGE:
# ./bundle_large_lifted_schur <sparse reconstruction file> <parameter file>
# params.txt : MAX_GMRES_ITERATIONS RESTARTS_FOR_GMRES TOLERANCE
./bundle_large_lifted_schur ~/Dataset/problem-138-19878-pre.txt ~/BundleAdjustment/codes/SSBA-4.0/params.txt
#./bundle_large_lifted_schur ~/Dataset/problem-49-7776-pre.txt
#./bundle_large ~/Dataset/problem-49-7776-pre.txt
