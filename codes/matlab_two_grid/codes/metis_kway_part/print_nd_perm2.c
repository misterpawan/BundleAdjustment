
#include "print_part.h"

int main(int argc, char *argv[]){

  idx_t      *part;
  MATCOO     ACOO;
  int nparts, err;
  FILE *pfile;
  FILE *f=NULL;
  
  f = fopen( argv[1], "r" );

  pfile = fopen("parts.txt", "w");

  
  nparts = atoi(argv[2]);
 // nparts = 4;

  if( f == NULL )
    {
      puts ( "cannot open file" ) ;
      exit(1) ;
    }

  err = readMatrixCoo(f, &ACOO);
 
// printf("nrows = %d\n", ACOO.M); 

  part = (idx_t*) calloc (ACOO.M, sizeof(idx_t));

  Interface_Metis_KWay_Part(nparts, &ACOO, part);

   for (int i=0; i<ACOO.M; i++)
    fprintf(pfile, "%d\n", (int)part[i]);

  //  fclose(pfile); fclose(f);
  // free(part);
  
  return 0;

}
