
clear
#before compiling for the first time
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/lib/intel64_lin:/home/shrutimoy.das/SuiteSparse/lib
#sudo ldconfig


#%% For UMFPACK
g++ -o mini_schur_complement mini_schur_complement.cpp -DMKL_LP64 -m64 -I ../SuiteSparse/include/ -I /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/ -L ../SuiteSparse/lib  -L /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin -lumfpack -lcxsparse -w -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#g++ -o mini_schur_complement mini_schur_complement.cpp -I ../SuiteSparse/include -L ../SuiteSparse/lib -lumfpack -lcxsparse

time ./mini_schur_complement

#%%For CHOLMOD
#g++ -o mini_schur_complement_CHOLMOD mini_schur_complement_CHOLMOD.cpp -I /usr/include/suitesparse -L /usr/lib/x86_64-linux-gnu -lcholmod -lumfpack -lspqr

#http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dcsrilu0_exampl1.c

#time ./mini_schur_complement_CHOLMOD


#test
#gcc -o dfgmres dfgmres.c -DMKL_LP64 -m64 -I /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/  -L /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin -w -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# for dfgmres test
#g++ -g -o dfgmres_test dfgmres_test.cpp -DMKL_LP64 -m64 -I ../SuiteSparse/include/ -I /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/ -L ../SuiteSparse/lib  -L /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin -lumfpack -lcxsparse -w -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl