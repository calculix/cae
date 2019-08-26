/*  testNDperm.c  */

#include "../misc.h"
#include "../../Graph.h"
#include "../../Perm.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------
   generate a nested dissection permutation object 
   for a n1 x n2 x n3 regular grid.

   created -- 96feb01, cca
   -----------------------------------------------
*/
{
Perm        *perm ;
FILE        *msgFile ;
int         j, msglvl, n1, n2, n3, nvtx ;
int         *newToOld, *oldToNew, *temp ;

if ( argc != 7 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile n1 n2 n3 permFile"
"\n        generate separators"
"\n msglvl   -- message level"
"\n msgFile  -- message file"
"\n n1       -- number of points in first direction"
"\n n2       -- number of points in second direction"
"\n n3       -- number of points in third direction"
"\n permFile -- file to contain the Perm structure"
"\n", argv[0]) ;
   return(1) ;
}
msglvl  = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n error in test4, file %s, line %d"
           "\n unable to open file %s", __FILE__, __LINE__, argv[2]) ;
   exit(-1) ; 
}
n1 = atoi(argv[3]) ;
n2 = atoi(argv[4]) ;
n3 = atoi(argv[5]) ;
/*
   -----------------------------
   create the permutation object
   -----------------------------
*/
nvtx = n1 * n2 * n3 ;
perm = Perm_new() ;
Perm_initWithTypeAndSize(perm, 3, nvtx) ;
newToOld = perm->newToOld ;
oldToNew = perm->oldToNew ;
/*
   ---------------------
   call the nd procedure
   ---------------------
*/
mkNDperm(n1, n2, n3, newToOld, 0, n1-1, 0, n2-1, 0, n3-1) ;
/*
   -----------------------------------------
   fill in the new-to-old permutation vector
   -----------------------------------------
*/
for ( j = 0 ; j < nvtx ; j++ ) {
   perm->oldToNew[perm->newToOld[j]] = j ;
}
if ( msglvl > 2 ) {
   Perm_writeForHumanEye(perm, msgFile) ;
}
/*
   ----------------------------------------
   write the permutation vector to its file
   ----------------------------------------
*/
if ( strcmp("none", argv[6]) != 0 ) {
   Perm_writeToFile(perm, argv[6]) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Perm_free(perm) ;

return(1) ; }

/*--------------------------------------------------------------------*/
