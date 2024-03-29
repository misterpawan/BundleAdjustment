/* Program to read a matrix, currently in IJV or triplet
or coordinate format, including matrix market file 

*** PARAMETERS ***

 f(INPUT) : file pointer to the matrix file
 ACOO(OUTPUT) : matrix in coordinate storage format (defined in mydefs.h)

 NOTE : Currently for matrix market files only, however, it may be 
 extended a bit for other formats like HB which are stored in compressed column
 or row formats

 *** NOTES ***

 Adapted from mmio.c from NIST by Pawan Kumar, Nov 24, 2009
 See http://math.nist.gov/MatrixMarket for details
 An extensive error handling is done by mmio.c while 
 reading or writing the matrix file. Most of the error handling part 
 is left on the mmio.c itself. Error hanling variable is ret_code
 
*/

//#include<stdio.h>
//#include<stdlib.h>
//#include "mydefs.h" 

//#include "newdefs.h"
//#include "brahm_defs.h"

#include "print_part.h"
#include "mmio.h"



//#include "mmio.h" /* Matrix market definitions (copied from NIST website !) */

int readMatrixCoo( FILE *f, MATCOO *ACOO )
{
  MM_typecode matcode; /* defined in mmio.h */
  //  int ret_code; /* error handler in mmio.c */
  int i;
  //  int err;
  //int nrows,ncols,nnz;
  int nargs;
  
  /* Read the banner */
  if(( mm_read_banner(f,&matcode))!=0){
    //    return -1;
  }
  
  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */
  
  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));exit(1);
  }
  
  /* find out size of sparse matrix .... */
  
  //  if ((ret_code = mm_read_mtx_crd_size(f, &ACOO->M, &ACOO->N, &ACOO->nnz)) !=0){
    // return ret_code; 
  //  printf("\n Error in reading size in readMatrix.c \n"); /* return an error code, if failed to read the size */
  //  exit(1);
  // }

  nargs = fscanf(f, "%d %d %d ", &ACOO->M, &ACOO->N, &ACOO->nnz );
  if(nargs != 3)
    {
      printf(" Could not read 3 arguments"); 
      exit(1);
    }

  
  /* reseve memory for matrices */

  ACOO->I = (int *) malloc(ACOO->nnz * sizeof(int)) ;
  ACOO->J = (int *) malloc(ACOO->nnz * sizeof(int)) ;
  ACOO->V = (double *) malloc(ACOO->nnz * sizeof(double));

  for (i=0; i<ACOO->nnz; i++)
    {
      nargs = fscanf(f, "%d %d %lg\n", &ACOO->I[i],
		     &ACOO->J[i], &ACOO->V[i]); /* Use lg or lf for double in
						// ANSI C */
      
      // f >> ACOO->I[i] >> ACOO->J[i] >> ACOO->V[i] ;
      if(nargs != 3)
	{
	  printf(" Could not read 3 arguments"); 
	  exit(1);
	}
     
        ACOO->I[i]--;ACOO->J[i]--;/* adjust from 1-based to 0-based */
  }
  
  if (f !=stdin) fclose(f);
  
  /* Write out matrix */
  //  mm_write_banner(stdout, matcode); 
//  mm_write_mtx_crd_size(stdout, ACOO->M, ACOO->N, ACOO->nnz);

  //  free(ACOO->I);  free(ACOO->J);  free(ACOO->V);

  return 0;

}
