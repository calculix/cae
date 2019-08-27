/*  scale.c  */

#include "../DenseMtx.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to scale a dense matrix object by a scalar
     A := alpha * A ;

   return values ---
      1 -- normal return
     -1 -- A is NULL
     -2 -- A has invalid type
     -3 -- alpha is NULL

   created -- 98nov06, cca
   -----------------------------------------------------
*/
int
DenseMtx_scale (
   DenseMtx   *A,
   double     alpha[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( A == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_scale()"
           "\n A is NULL\n") ;
   return(-1) ;
}
if ( A->type != SPOOLES_REAL && A->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n error in DenseMtx_scale()"
           "\n A has invalid type\n") ;
   return(-2) ;
}

if ( alpha == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_scale()"
           "\n alpha is NULL\n") ;
   return(-3) ;
}
if ( A->type == SPOOLES_REAL ) {
   double   ralpha = alpha[0] ;
   if ( ralpha == 1.0 ) {
      return(1) ;
   } else {
      double   *entries ;
      int      colinc, ncol, nrow, rowinc ;

      entries = DenseMtx_entries(A) ;
      rowinc  = DenseMtx_rowIncrement(A) ;
      colinc  = DenseMtx_columnIncrement(A) ;
      DenseMtx_dimensions(A, &nrow, &ncol) ;
      if (  (rowinc == 1 && colinc == nrow) 
         || (rowinc == ncol && colinc == 1) ) {
         if ( ralpha == 0.0 ) {
            DVzero(nrow*ncol, entries) ;
         } else {
            DVscale(nrow*ncol, entries, ralpha) ;
         }
      } else {
         int   ii, jj ;
         if ( ralpha == 0.0 ) {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               for ( ii = 0 ; ii < nrow ; ii++ ) {
                  entries[ii*rowinc + jj*colinc] = 0.0 ;
               }
            }
         } else {
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               for ( ii = 0 ; ii < nrow ; ii++ ) {
                  entries[ii*rowinc + jj*colinc] *= ralpha ;
               }
            }
         }
      }
   }
} else if ( A->type == SPOOLES_COMPLEX ) {
   double ralpha = alpha[0], ialpha = alpha[1] ;
   if ( ralpha == 1.0 && ialpha == 0.0 ) {
      return(1) ;
   } else {
      double   *entries ;
      int      colinc, ncol, nrow, rowinc ;

      entries = DenseMtx_entries(A) ;
      rowinc  = DenseMtx_rowIncrement(A) ;
      colinc  = DenseMtx_columnIncrement(A) ;
      DenseMtx_dimensions(A, &nrow, &ncol) ;
      if (  (rowinc == 1 && colinc == nrow) 
         || (rowinc == ncol && colinc == 1) ) {
         if ( ralpha == 0.0 && ialpha == 0.0 ) {
            ZVzero(nrow*ncol, entries) ;
         } else {
            ZVscale(nrow*ncol, entries, ralpha, ialpha) ;
         }
      } else {
         if ( ralpha == 0.0 && ialpha == 0.0 ) {
            int      ii, jj, off ;
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               for ( ii = 0 ; ii < nrow ; ii++ ) {
                  off = ii*rowinc + jj*colinc ;
                  entries[2*off]   = 0.0 ;
                  entries[2*off+1] = 0.0 ;
               }
            }
         } else {
            double   yi, yr ;
            int      ii, jj, off ;
            for ( jj = 0 ; jj < ncol ; jj++ ) {
               for ( ii = 0 ; ii < nrow ; ii++ ) {
                  off = ii*rowinc + jj*colinc ;
                  yr = entries[2*off] ; yi = entries[2*off+1] ;
                  entries[2*off]   = yr * ralpha - yi * ialpha ;
                  entries[2*off+1] = yr * ialpha - yi * ralpha ;
               }
            }
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
