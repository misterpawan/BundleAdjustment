
clear

#before compiling for the first time
#LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64_lin
#sudo ldconfig


#%% For UMFPACK
g++ -o mini_schur_complement mini_schur_complement.cpp -DMKL_ILP64 -m64 -I /usr/include/suitesparse -I /opt/intel/mkl/include -L /usr/lib/x86_64-linux-gnu -L /opt/intel/mkl/lib/intel64_lin -lumfpack -lcxsparse -W -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

time ./mini_schur_complement

#%%For CHOLMOD
#g++ -o mini_schur_complement_CHOLMOD mini_schur_complement_CHOLMOD.cpp -I /usr/include/suitesparse -L /usr/lib/x86_64-linux-gnu -lcholmod -lumfpack -lspqr

#http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dcsrilu0_exampl1.c

#time ./mini_schur_complement_CHOLMOD