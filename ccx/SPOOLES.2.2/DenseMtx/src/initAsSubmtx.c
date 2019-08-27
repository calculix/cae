/*  initAsSubmtx.c  */

#include "../DenseMtx.h"

#define MYDEBUG 1

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- initialize as a submatrix of another DenseMtx object.
         B = A(firstrow:lastrow, firstcol:lastcol)
      note, B only points into the storage of A.

   return values --
      1 -- normal return
     -1 -- B is NULL
     -2 -- A is NULL
     -3 -- A has invalid type
     -4 -- requested rows are invalid
     -5 -- requested columns are invalid

   created -- 98nov11, cca
   ----------------------------------------------------------------
*/
int
DenseMtx_initAsSubmatrix (
   DenseMtx   *B,
   DenseMtx   *A,
   int        firstrow,
   int        lastrow,
   int        firstcol,
   int        lastcol
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  B == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_initAsSubmatrix()"
           "\n B is NULL\n") ;
   return(-1) ;
}
if (  A == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_initAsSubmatrix()"
           "\n A is NULL\n") ;
   return(-2) ;
}
if ( A->type != SPOOLES_REAL && A->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n error in DenseMtx_initAsSubmatrix()"
           "\n invalid type %d\n", A->type) ;
   return(-3) ;
}
if (  firstrow < 0 || lastrow >= A->nrow ) {
   fprintf(stderr, "\n error in DenseMtx_initAsSubmatrix()"
           "\n %d rows in A, firstrow is %d\n", A->nrow, firstrow) ;
   return(-4) ;
}
if (  firstcol < 0 || lastcol >= A->ncol ) {
   fprintf(stderr, "\n error in DenseMtx_initAsSubmatrix()"
           "\n %d columns in A, firstcol is %d\n", A->ncol, firstcol) ;
   return(-5) ;
}
/*
   ---------------------
   set the scalar fields
   ---------------------
*/
B->type    = A->type ;
B->rowid   = A->rowid ;
B->colid   = A->colid ;
B->nrow    = lastrow - firstrow + 1 ;
B->ncol    = lastcol - firstcol + 1 ;
B->inc1    = A->inc1  ;
B->inc2    = A->inc2  ;
B->rowind  = A->rowind + firstrow ;
B->colind  = A->colind + firstcol ;
if ( A->type == SPOOLES_REAL ) {
   B->entries = A->entries + firstrow*A->inc1 + firstcol*A->inc2 ;
} else {
   B->entries = A->entries + 2*(firstrow*A->inc1 + firstcol*A->inc2) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
