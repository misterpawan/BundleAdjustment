// function to perform x = x + beta*A(:,j)
// adapted from CSparse of T. Davis of UoF
// Pawan kumar, Inria futurs, Jan 3, 2010

// #include "brahm_defs.h"
#include "print_part.h"

int Brahm_scatter(MATCSC *A, int j, double beta, 
		  int *w, double *x, int mark,
		  MATCSC *C, int nz){
  
  int i, p, *Ap, *Ai, *Ci ;
  double *Ax ;
  
  Ap = A->J ; Ai = A->I ; Ax = A->V ; Ci = C->I ;
  
  for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      i = Ai [p] ;                            /* A(i,j) is nonzero */
      if (w [i] < mark)
        {
	  w [i] = mark ;                      /* i is new entry in column j */
	  Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
	  if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
      else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
  return (nz) ;
}
