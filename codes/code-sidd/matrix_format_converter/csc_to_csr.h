/*

Conversion from Compressed Sparse Column to Compressed Sparse Row format 

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


struct csr_matrix_t
{
  
  int m;   // Number of rows in the matrix 

  int n;   // Number of columns in the matrix 

  int nnz; // Number of stored (nonzero) entries

  void* values;   // Array of stored (nonzero) entries of the matrix 

  int* colidx;  // Array of column indices of the stored (nonzero) entries of the matrix 

  int* rowptr;  // Array of indices into the colidx and values arrays, for each column

};




struct csr_matrix_t*
csc_to_csr (struct csc_matrix_t* A)
{
  int *row_nnz = NULL; /* # of nonzeros in each row */
  int i, j;
  struct csr_matrix_t* B = calloc (1, sizeof (struct csr_matrix_t));

  printf("--- start csc to csr conversion ---\n");
  
  B->m = A->m;              // Number of Rows
  B->n = A->n;              // Number of Columns
  B->nnz = A->nnz;          // Number of Non-zeros

  B->rowptr = malloc ((A->m + 1) * sizeof (int));
  B->colidx = malloc (A->nnz * sizeof (int));
  B->values = malloc (A->nnz * sizeof (double));
  

  /* count # of non-zeros per row */

  row_nnz = calloc (A->m, sizeof(int));

  for (i = 0; i < A->nnz; i++)
    {
      int k = A->rowidx[i];
      row_nnz[k]++;
    }

  /*
   *  initialize CSR's row pointers.
   *  reset row_nnz to zero.
   */

  B->rowptr[0] = 0;
  for (i = 1; i <= A->m; i++)
    {
      B->rowptr[i] = B->rowptr[i-1] + row_nnz[i-1];
      row_nnz[i-1] = 0;
    }

  
  /*
   *  convert from CSC to CSR.
   *  use row_nnz to keep track of the number of non-zeros
   *  added to each row.
   */
    
  
  double* in_values  = (double*) (A->values);
  double* out_values = (double*) (B->values);

  for (j = 0; j < A->n; j++)
	{
	   int k;
	   int nnz_col;    /* # of non-zeros in column j of A */

	   nnz_col = A->colptr[j+1] - A->colptr[j];

	   for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	   {
	       int i = A->rowidx[ k ];             /* row index */
	       double a = in_values[ k ];          /* non-zero value */
	       int h = B->rowptr[i] + row_nnz[i];  /* non-zero position */

	       /* add the non-zero A(i,j) to B */

	       B->colidx[ h ] = j;
	       out_values[ h ] = a;

	       row_nnz[i]++;
	   }
	}
    

  free (row_nnz);    

  printf("--- finished csc to csr conversion  ---\n");

  return B;

}
