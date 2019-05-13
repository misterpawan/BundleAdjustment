# your FORTRAN compiler
F90=gfortran

# your C compiler
C=gcc

# set here link reference needed for BLAS & LAPACK
BLASLAPACK=-L/usr/lib -llapack -lblas

# compilation options
opt = -Ofast #-fdefault-integer-8

# directory where AGMG sources are located
agmgdir=../SRC

lib=$(BLASLAPACK)
libc=-lgfortran -lm $(BLASLAPACK)

export F90
export opt

list=Example_seq.o $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o
list_c=Example_seq_c.o $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o

default: Example_seq

all: Example_seq Example_seq_c

Example_seq: Example_seq.o
	cd $(agmgdir);make dseq
	$(F90) $(opt) -o Example_seq $(list) $(lib)

Example_seq.o: Example_seq.f90
	$(F90) $(opt) -c Example_seq.f90

#Example_seq_c: Example_seq_c.o
#	cd $(agmgdir);make dseq
#	$(F90) $(opt) -o Example_seq_c $(list_c) $(lib)

Example_seq_c: Example_seq_c.o
	cd $(agmgdir);make dseq
	$(C) $(opt) -o Example_seq_c $(list_c) $(libc)

Example_seq_c.o: Example_seq_c.c
	$(C) $(opt) -c Example_seq_c.c
