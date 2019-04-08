dseq:dagmg.o dagmg_mumps.o
sseq:sagmg.o sagmg_mumps.o
cseq:cagmg.o cagmg_mumps.o
zseq:zagmg.o zagmg_mumps.o
dmpi:dagmg_mpi.o
smpi:sagmg_mpi.o
cmpi:cagmg_mpi.o
zmpi:zagmg_mpi.o
lapack:lapack_agmg.o
blas:blas_agmg.o

dagmg.o: dagmg.f90
	$(F90) $(opt) -c dagmg.f90
dagmg_mumps.o: dagmg_mumps.f90
	$(F90) $(opt) -c dagmg_mumps.f90
dagmg_mpi.o: dagmg_mpi.f90
	$(F90) $(opt) $(incl) -c dagmg_mpi.f90

sagmg.o: sagmg.f90
	$(F90) $(opt) -c sagmg.f90
sagmg_mumps.o: sagmg_mumps.f90
	$(F90) $(opt) -c sagmg_mumps.f90
sagmg_mpi.o: sagmg_mpi.f90
	$(F90) $(opt) $(incl) -c sagmg_mpi.f90

cagmg.o: cagmg.f90
	$(F90) $(opt) -c cagmg.f90
cagmg_mumps.o: cagmg_mumps.f90
	$(F90) $(opt) -c cagmg_mumps.f90
cagmg_mpi.o: cagmg_mpi.f90
	$(F90) $(opt) $(incl) -c cagmg_mpi.f90

zagmg.o: zagmg.f90
	$(F90) $(opt) -c zagmg.f90
zagmg_mumps.o: zagmg_mumps.f90
	$(F90) $(opt) -c zagmg_mumps.f90
zagmg_mpi.o: zagmg_mpi.f90
	$(F90) $(opt) $(incl) -c zagmg_mpi.f90

lapack_agmg.o: lapack_agmg.f
	$(F90) $(opt) -c lapack_agmg.f
blas_agmg.o: blas_agmg.f
	$(F90) $(opt) -c blas_agmg.f
