/*  permute.c  */

#include "../DenseMtx.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to permute the rows of an object

   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_permuteRows (
   DenseMtx   *mtx,
   IV         *oldToNewIV
) {
A2    a2 ;
int   ii, irow, maxnrow, nrow ;
int   *oldToNew, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || oldToNewIV == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_permuteRows(%p,%p)"
           "\n bad input\n", mtx, oldToNewIV) ;
   exit(-1) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( nrow <= 0 ) {
   return ;
}
/*
   ----------------------------------------------
   overwrite the old row ids with the new row ids
   ----------------------------------------------
*/
IV_sizeAndEntries(oldToNewIV, &maxnrow, &oldToNew) ;
for ( ii = 0 ; ii < nrow ; ii++ ) {
   irow = rowind[ii] ;
   if ( irow < 0 || irow >= maxnrow ) {
      fprintf(stderr, "\n fatal error in DenseMtx_permuteRows(%p,%p)"
              "\n irow = %d, maxnrow = %d", 
              mtx, oldToNewIV, irow, maxnrow) ;
      exit(-1) ;
   }
   rowind[ii] = oldToNew[rowind[ii]] ;
}
/*
   ------------------------------------
   now sort the rows in ascending order
   ------------------------------------
*/
A2_setDefaultFields(&a2) ;
DenseMtx_setA2(mtx, &a2) ;
A2_sortRowsUp(&a2, nrow, rowind) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to permute the columns of an object

   created -- 98may02, cca
   ----------------------------------------------
*/
void
DenseMtx_permuteColumns (
   DenseMtx   *mtx,
   IV          *oldToNewIV
) {
A2    a2 ;
int   ii, jcol, maxncol, ncol ;
int   *oldToNew, *colind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || oldToNewIV == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_permuteColumns(%p,%p)"
           "\n bad input\n", mtx, oldToNewIV) ;
   exit(-1) ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
if ( ncol <= 0 ) {
   return ;
}
/*
   ----------------------------------------------------
   overwrite the old column ids with the new column ids
   ----------------------------------------------------
*/
IV_sizeAndEntries(oldToNewIV, &maxncol, &oldToNew) ;
for ( ii = 0 ; ii < ncol ; ii++ ) {
   jcol = colind[ii] ;
   if ( jcol < 0 || jcol >= maxncol ) {
      fprintf(stderr, 
              "\n fatal error in DenseMtx_permuteColumns(%p,%p)"
              "\n jcol = %d, maxncol = %d", 
              mtx, oldToNewIV, jcol, maxncol) ;
      exit(-1) ;
   }
   colind[ii] = oldToNew[jcol] ;
}
/*
   ------------------------------------
   now sort the rows in ascending order
   ------------------------------------
*/
A2_setDefaultFields(&a2) ;
DenseMtx_setA2(mtx, &a2) ;
A2_sortColumnsUp(&a2, ncol, colind) ;

return ; }

/*--------------------------------------------------------------------*/
