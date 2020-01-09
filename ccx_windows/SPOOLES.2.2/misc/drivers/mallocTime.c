#include <stdio.h>
#include "../../IV.h"
#include "../../Utilities.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   this program tries to get some idea of the malloc()/free()
   characteristics of a machine and operating system. it was
   written to try and understand the initialization times for
   the FrontMtx object for the i4a matrix with 99192 rows and
   columns.

   the program reads in a vector that contains the byte counts
   for the factor submatrices. depending on the input parameters,
   the program allocates the storage and then free'ing the storage.

   zeroflag -- 
          0     --> do not zero the storage
      otherwise --> zero the storage
   sortflag -- 
          0     --> do not sort the byte counts in ascending order
      otherwise --> sort the byte counts in ascending order

   created -- 98sep05, cca
   ---------------------------------------------------------------
*/
int
main ( int argc, char *argv[] ) {
double   t1, t2 ;
double   *dvec ;
double   **pdvecs ;
int      item, nitem, sortflag, sum, zeroflag ;
int      *nbytes ;
IV       *nbytesIV ;

if ( argc != 3 ) {
   fprintf(stdout, "\n usage : a.out zeroflag sortflag ") ;
   return(1) ;
}
zeroflag = atoi(argv[1]) ;
sortflag = atoi(argv[2]) ;
if ( zeroflag == 0 ) {
   fprintf(stdout, "\n storage not zero'd") ;
} else {
   fprintf(stdout, "\n storage zero'd") ;
}
if ( sortflag == 0 ) {
   fprintf(stdout, "\n byte counts not sorted") ;
} else {
   fprintf(stdout, "\n byte counts sorted") ;
}
/*
   -------------------------------------------------
   read in the vector that contains the bytes counts 
   for the SubMtx objects from the i4a matrix
   -------------------------------------------------
*/
nbytesIV = IV_new() ;
IV_readFromFile(nbytesIV, "nbytes.ivf") ;
IV_sizeAndEntries(nbytesIV, &nitem, &nbytes) ;
sum = IV_sum(nbytesIV) ;
fprintf(stdout, "\n %d items read in, sum = %d",
        nitem, sum) ;
fflush(stdout) ;
if ( sortflag != 0 ) {
   IVqsortUp(nitem, nbytes) ; 
}
/*
   ------------------------------------
   now time the malloc()'s and free()'s
   ------------------------------------
*/
pdvecs = PDVinit(nitem) ;
if ( zeroflag == 0 ) {
   MARKTIME(t1) ;
   for ( item = 0 ; item < nitem ; item++ ) {
      pdvecs[item] = DVinit2(nbytes[item]/sizeof(double)) ;
   }
   MARKTIME(t2) ;
   fprintf(stdout, "\n  CPU %10.5f for DVinit2", t2 - t1) ;
   
   MARKTIME(t1) ;
   for ( item = 0 ; item < nitem ; item++ ) {
      DVfree(pdvecs[item]) ;
   }
   MARKTIME(t2) ;
   fprintf(stdout, "\n  CPU %10.5f for DVfree", t2 - t1) ;
} else {
   MARKTIME(t1) ;
   for ( item = 0 ; item < nitem ; item++ ) {
      pdvecs[item] = DVinit(nbytes[item]/sizeof(double), 0.0) ;
   }
   MARKTIME(t2) ;
   fprintf(stdout, "\n  CPU %10.5f for DVinit", t2 - t1) ;
   
   MARKTIME(t1) ;
   for ( item = 0 ; item < nitem ; item++ ) {
      DVfree(pdvecs[item]) ;
   }
   MARKTIME(t2) ;
   fprintf(stdout, "\n  CPU %10.5f for DVfree", t2 - t1) ;
}
fprintf(stdout, "\n") ;
fflush(stdout) ;
return(1) ; }

/*--------------------------------------------------------------------*/
