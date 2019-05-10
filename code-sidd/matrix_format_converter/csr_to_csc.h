/*

Conversion from Compressed Sparse Row to Compressed Sparse Column format 

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


struct csc_matrix_t
{
  
  int m;  // Number of rows in the matrix 

  int n;  // Number of columns in the matrix
  
  int nnz;  // Number of stored (nonzero) entries

  void* values;  // Array of stored (nonzero) entries of the matrix

  int* rowidx; // Array of row indices of the stored (nonzero) entries of the matrix 

  int* colptr;  // Array of indices into the rowidx and values arrays, for each column 

};




struct csc_matrix_t*
csr_to_csc (const struct csr_matrix_t* A)
{
  int *col_nnz = NULL;   // # of nonzeros in each column 
  int i;
  struct csc_matrix_t* B = calloc (1, sizeof (struct csc_matrix_t));

  printf("--- start csr to csc conversion ---\n");

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;

  B->colptr = malloc ((A->n + 1) * sizeof (int));
  B->rowidx = malloc (A->nnz * sizeof (int));
  B->values = malloc (A->nnz * sizeof (double));


  /* count # of non-zeros per column */

  col_nnz = calloc (A->n, sizeof(int));

  for (i = 0; i < A->nnz; i++)
   {
      int k = A->colidx[i];
      col_nnz[k]++;
   }

  /*
   *  initialize CSC's column pointers.
   *  reset col_nnz to zero.
   */

  B->colptr[0] = 0;
  for (i = 1; i <= A->n; i++)
   {
      B->colptr[i] = B->colptr[i-1] + col_nnz[i-1];
      col_nnz[i-1] = 0;
   }

  /*
   *  convert from CSR to CSC.
   *  use col_nnz to keep track of the number of non-zeros
   *  added to each column.
   */
  
    
  double* in_values  = (double*) (A->values);
  double* out_values = (double*) (B->values);

  for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    // # of non-zeros in row i of A 

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];             // column index 
	      double a = in_values[ k ];          // non-zero value 
	      int h = B->colptr[j] + col_nnz[j];  // non-zero position 

	      /* add the non-zero A(i,j) to B */

	      B->rowidx[ h ] = i;
	      out_values[ h ] = a;

	      col_nnz[j]++;
	    }
	}
    
  

  free (col_nnz);  

  printf("--- finished csr to csc conversion  ---\n");

  return B;

}
