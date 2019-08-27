/*  makeStaircase.c  */

#include "../A2.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   purpose -- to permute the rows so the matrix is in staircase form

   created -- 98may25, cca
   -----------------------------------------------------------------
*/
void
A2_makeStaircase (
   A2   *mtxA
) {
double   imag, real, value ;
int      irow, jcol, ncol, nrow ;
int      *firstnonzero ;
/*
   --------------
   check the data
   --------------
*/
if ( mtxA == NULL ) {
   fprintf(stderr, "\n fatal error in A2_staircase(%p)"
           "\n bad input\n", mtxA) ;
   exit(-1) ;
}
nrow = A2_nrow(mtxA) ;
ncol = A2_ncol(mtxA) ;
/*
   ---------------------------------------------
   fill firstnonzero[irow] with the first column 
   that contains a nonzero entry in this row.
   ---------------------------------------------
*/
firstnonzero = IVinit(nrow, -1) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      if ( A2_IS_REAL(mtxA) ) {
         A2_realEntry(mtxA, irow, jcol, &value) ;
         if ( value != 0.0 ) {
            break ;
         }
      } else if ( A2_IS_COMPLEX(mtxA) ) {
         A2_complexEntry(mtxA, irow, jcol, &real, &imag) ;
         if ( real != 0.0 || imag != 0.0 ) {
            break ;
         }
      }
   }
   firstnonzero[irow] = jcol ;
}
/*
   ---------------------------------------------------
   sort the rows in the order of their leading nonzero
   ---------------------------------------------------
*/
A2_sortRowsUp(mtxA, nrow, firstnonzero) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(firstnonzero) ;

return ; }

/*--------------------------------------------------------------------*/
