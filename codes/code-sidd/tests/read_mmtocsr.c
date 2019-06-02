/* reads matrix market format (coordinate) and returns
   csr format */

#include "my_util.h"
#include "mmio.h"
#include "matrix_io.h"

void read_mm_matrix (char *fn, int *m, int *n, int *nz,
      int **i_idx, int **j_idx, double **a);
void coo2csr_in(int n, int nz, double *a, int *i_idx, int *j_idx);


int main(int argc, char *argv[])
{

  f = fopen(argv[1], "r");     // Read the matrix during compile time
  MM_typecode matcode;
  FILE *f;
  int i,k;
  int base;

  if ((f = fopen(fn, "r")) == NULL) {
    printf ("can't open file <%s> \n", fn);
    exit(1);
  }
  if (mm_read_banner(f, &matcode) != 0){
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (! (mm_is_matrix(matcode) && mm_is_sparse(matcode)) ){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    exit(1);
  }

  /* find out size of sparse matrix .... */

  fscanf(f, "%*s %d %d %d", m, n, nz);    

  fscanf(f, "%*s %d", &base);     /* indentifies starting index */

  /* reserve memory for matrices */
  if (mm_is_symmetric(matcode)){
    *i_idx = (int *) my_malloc(*nz *2 * sizeof(int));
    *j_idx = (int *) my_malloc(*nz *2 * sizeof(int));
    *a = (double *) my_malloc(*nz *2 * sizeof(double));
  }
  else {
    *i_idx = (int *) my_malloc(*nz * sizeof(int));
    *j_idx = (int *) my_malloc(*nz * sizeof(int));
    *a = (double *) my_malloc(*nz * sizeof(double));
  }

  if (!(*i_idx) || !(*j_idx) || !(*a)){
    printf ("cannot allocate memory for %d, %d, %d sparse matrix\n", *m, *n, *nz);
    exit(1);
  }


  k=0;
  for (i=0; i<*nz; i++)  {
    if (mm_is_pattern(matcode)){
      fscanf(f, "%d %d", &(*i_idx)[i], &(*j_idx)[i]);
      (*i_idx)[i] -= base;  /* adjust from 1-based to 0-based */
      (*j_idx)[i] -= base;

      (*a)[i] = random_double(-1, 1);
    }
    else if (mm_is_real(matcode)){
      fscanf(f, "%d %d %lg", &(*i_idx)[i], &(*j_idx)[i], &(*a)[i]);
      (*i_idx)[i] -= base;  /* adjust from 1-based to 0-based */
      (*j_idx)[i] -= base;
    }

    if (mm_is_symmetric(matcode)){
      if ( (*i_idx)[i] != (*j_idx)[i] ){
  (*i_idx)[*nz+k] = (*j_idx)[i];
  (*j_idx)[*nz+k] = (*i_idx)[i];
  (*a)[*nz+k] = (*a)[i];
  k++;
      }
    }
  }
  *nz += k;

  fclose(f);

  coo2csr_in (*m, *nz, *a, *i_idx, *j_idx);
}


/* converts COO format to CSR format, in-place,
   if SORT_IN_ROW is defined, each row is sorted in column index.
   On return, i_idx contains row_start position */

void coo2csr_in(int n, int nz, double *a, int *i_idx, int *j_idx)
{
  int *row_start;
  int i, j;
  int init, i_next, j_next, i_pos;
  double dt, a_next;

  row_start = (int *)malloc((n+1)*sizeof(int));
  if (!row_start){
    printf ("coo2csr_in: cannot allocate temporary memory\n");
    exit (1);
  }
  for (i=0; i<=n; i++) row_start[i] = 0;

  /* determine row lengths */
  for (i=0; i<nz; i++) row_start[i_idx[i]+1]++;

  for (i=0; i<n; i++) row_start[i+1] += row_start[i];

  for (init=0; init<nz; ){
    dt = a[init];
    i = i_idx[init];
    j = j_idx[init];
    i_idx[init] = -1;
    while (1){
      i_pos = row_start[i];
      a_next = a[i_pos];
      i_next = i_idx[i_pos];
      j_next = j_idx[i_pos];

      a[i_pos] = dt;
      j_idx[i_pos] = j;
      i_idx[i_pos] = -1;
      row_start[i]++;
      if (i_next < 0) break;
      dt = a_next;
      i = i_next;
      j = j_next;

    }
    init++;
    while ((i_idx[init] < 0) && (init < nz))  init++;
  }


  /* shift back row_start */
  for (i=0; i<n; i++) i_idx[i+1] = row_start[i];
  i_idx[0] = 0;


  for (i=0; i<n; i++){
    sort (j_idx, a, i_idx[i], i_idx[i+1]);
  }

}
