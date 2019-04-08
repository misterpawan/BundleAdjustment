% make_agmg_octfile:
%  script that generates agmg.oct from source files;
%  agmg.oct implements the agmg function; enter 'help agmg' for usage.
%  The file 'agmg.oct' should be moved to a directory in Octave PATH.
%  Once agmg.oct is built, all other files from the AGMG package may 
%  be deleted.
%
%  If you have also the file 'agmg.m' (needed for the agmg function
%  withing Matlab), make sure that agmg.oct is either in the same 
%  directory, or in a directory that have precedence in Octave PATH.
%
mkoctfile agmg.cc dagmg_oct.f90 zagmg_oct.f90 dagmg_mumps.f90 zagmg_mumps_nocommonpart.f90 -lgfortran -lm -llapack -lblas
