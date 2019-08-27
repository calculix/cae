/*  test1.c  */

#include <stdio.h>

/*--------------------------------------------------------------------*/

int
main (
   int    argc,
   char   *argv[]
) {
char   ctemp ;
FILE   *fp ;
int    itemp, rc ;

if ( argc != 2 ) {
   fprintf(stdout, "\n usage : filename") ;
   return(0) ;
}
if ( (fp = fopen(argv[1], "r")) == NULL ) {
   fprintf(stderr, "\n unable to open file %s\n", argv[1]) ;
}
/*
   ---------------
   check the input
   ---------------
*/
while ( 1 ) {
   rc = fscanf(fp, "%d%c", &itemp, &ctemp) ;
   if ( rc != 2 ) {
      fprintf(stdout, "\n error, rc = %d", rc) ;
      break ;
   }
   fprintf(stdout, "\n itemp = '%d', ctemp = '%c'", itemp, ctemp) ;
   if ( ctemp == '\n' ) {
      fprintf(stdout, "\n newline") ;
   }
   if ( ctemp == EOF ) {
      fprintf(stdout, "\n end of file") ;
      break ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
