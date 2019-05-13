/* *******************************************************************/
/*   Convert matrix in coordinate storage format to compressed column
     storage format */

/*******************************************************************/
/*
   This function converts the coordinate storage format to 
 compressed column format. Adapted from SparseLib++ (Free
 Licence) by Pawan Kumar on Nov 24, 2009 

 *** PARAMETERS ***

 ACOO (INPUT) : Input matrix in coordinate storage format
 ACSC (OUTPUT) : Output matrix in CSC format
*/

//#include<stdio.h>
//#include "mydefs.h" /* MATCOO and MATCSC defined here */

// #include "brahm_defs.h"

#include "print_part.h"


int coocsc3( int ncol, int nnz, double *val, int *rowind, int *colind, double *val_new, int *rowind_new, int *colptr)
{
  int i,j;
  int *tally;

  //   printf("\n ncol : %d, nnz : %d \n", ncol, nnz );
  tally = (int*)malloc((ncol+1)*sizeof(int));
    
  for(i=0;i<ncol+1;i++) tally[i]=0;
  
  for (i=0;i<nnz;i++) tally[colind[i]]++;
  colptr[0] = 0;

  for (j=0;j<ncol;j++) colptr[j+1] = colptr[j]+tally[j];
  // Make copy of rowptr for use in second pass.
  for (i=0;i<ncol+1;i++) tally[i] = colptr[i];
  // Second pass through nonzeros. Fill in index and value entries.
  for (i=0;i<nnz;i++)
    {
      val_new[tally[colind[i]]] = val[i];
      rowind_new[tally[colind[i]]] = rowind[i];
      tally[colind[i]]++;
    }

  free(tally);
  return 0;
}

/*** end of function writeMatrixCsc ***/
