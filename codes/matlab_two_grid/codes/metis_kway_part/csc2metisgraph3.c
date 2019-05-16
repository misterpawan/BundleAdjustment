//#include<metis.h> // necessary to include : to include all
//#later!

// #include "brahm_defs.h"

#include "print_part.h"

int csc2metisgraph3 (int n, int nnz, int *ir, int *jc, double *pr, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt)
{
  int i, j, jbar;

  /* Scan the matrix, not copying diagonal terms, and rounding doubles
   * to integer weights */

  xadj[0] = 0;
  jbar = 0; 
  for (i = 1; i <= n; i++) 
    {
       for (j = jc[i-1]; j < jc[i]; j++) 
	{
	  if (ir[j] != i-1) {
	    adjncy[jbar] = ir[j];

	  //  printf("\n jbar= %d \n", jbar);

	    adjwgt[jbar] = (int) pr[j];
	    jbar++;
	  } else {
	    vwgt[i-1] = (int) pr[j];
	  }
	}
      xadj[i] = jbar;
    }
  return 0;

}