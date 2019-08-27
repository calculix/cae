/*  adjForATA.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return an IVL object that contains 
   the adjacency structure of A^TA.

   created -- 98jan28, cca
   ----------------------------------
*/
IVL *
InpMtx_adjForATA (
   InpMtx   *inpmtxA
) {
InpMtx   *inpmtxATA ;
int      firstcol, firstrow, irow, jvtx, lastcol, lastrow,
         loc, ncol, nent, nrow, size ;
int      *ind, *ivec1, *ivec2 ;
IVL      *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtxA == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_adjForATA(%p)"
           "\n NULL input\n", inpmtxA) ;
   exit(-1) ;
}
/*
   ----------------------------------------------------------
   change the coordinate type and storage mode to row vectors
   ----------------------------------------------------------
*/
InpMtx_changeCoordType(inpmtxA, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(inpmtxA, INPMTX_BY_VECTORS) ;
nent     = InpMtx_nent(inpmtxA) ;
ivec1    = InpMtx_ivec1(inpmtxA) ;
ivec2    = InpMtx_ivec2(inpmtxA) ;
firstrow = IVmin(nent, ivec1, &loc) ;
lastrow  = IVmax(nent, ivec1, &loc) ;
firstcol = IVmin(nent, ivec2, &loc) ;
lastcol  = IVmax(nent, ivec2, &loc) ;
if ( firstrow < 0 || firstcol < 0 ) {
   fprintf(stderr, "\n fatal error"
           "\n firstrow = %d, firstcol = %d"
           "\n lastrow  = %d, lastcol  = %d",
           firstrow, firstcol, lastrow, lastcol) ;
   exit(-1) ;
}
nrow = 1 + lastrow ;
ncol = 1 + lastcol ;
/*
   -----------------------------------------------------------
   create the new InpMtx object to hold the structure of A^TA
   -----------------------------------------------------------
*/
inpmtxATA = InpMtx_new() ;
InpMtx_init(inpmtxATA, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, 0, 0) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   InpMtx_vector(inpmtxA, irow, &size, &ind) ;
   InpMtx_inputMatrix(inpmtxATA, size, size, 1, size, ind, ind) ;
}
for ( jvtx = 0 ; jvtx < nrow ; jvtx++ ) {
   InpMtx_inputEntry(inpmtxATA, jvtx, jvtx) ;
}
InpMtx_changeStorageMode(inpmtxATA, INPMTX_BY_VECTORS) ;
/*
   -------------------
   fill the IVL object
   -------------------
*/
adjIVL = IVL_new() ;
IVL_init1(adjIVL, IVL_CHUNKED, nrow) ;
for ( jvtx = 0 ; jvtx < ncol ; jvtx++ ) {
   InpMtx_vector(inpmtxATA, jvtx, &size, &ind) ;
   IVL_setList(adjIVL, jvtx, size, ind) ;
}
/*
   ------------------------------
   free the working InpMtx object
   ------------------------------
*/
InpMtx_free(inpmtxATA) ;

return(adjIVL) ; }

/*--------------------------------------------------------------------*/
