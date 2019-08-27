/*  mmm.c  */

#include "../Pencil.h"

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   compute Y := Y + (A + sigma*B)*X

   created -- 98may02, cca
   --------------------------------
*/
void
Pencil_mmm (
   Pencil     *pencil,
   DenseMtx   *Y,
   DenseMtx   *X
) {
int   ncolX, ncolY, nrhs, nrow, nrowX, nrowY ;
/*
   ---------------
   check the input
   ---------------
*/
if ( pencil == NULL || Y == NULL || X == NULL ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n bad input\n", pencil, Y, X) ;
   exit(-1) ;
}
if ( !(PENCIL_IS_REAL(pencil) || PENCIL_IS_COMPLEX(pencil)) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n bad type %d for pencil\n", pencil, Y, X, pencil->type) ;
   exit(-1) ;
}
if ( !(DENSEMTX_IS_REAL(Y) || DENSEMTX_IS_COMPLEX(Y)) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n bad type %d for Y\n", pencil, Y, X, Y->type) ;
   exit(-1) ;
}
if ( !(DENSEMTX_IS_REAL(X) || DENSEMTX_IS_COMPLEX(X)) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n bad type %d for X\n", pencil, Y, X, X->type) ;
   exit(-1) ;
}
if ( PENCIL_IS_REAL(pencil) && !DENSEMTX_IS_REAL(Y) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n pencil is real, Y is not\n", pencil, Y, X) ;
   exit(-1) ;
}
if ( PENCIL_IS_REAL(pencil) && !DENSEMTX_IS_REAL(X) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n pencil is real, X is not\n", pencil, Y, X) ;
   exit(-1) ;
}
if ( PENCIL_IS_COMPLEX(pencil) && !DENSEMTX_IS_COMPLEX(Y) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n pencil is complex, Y is not\n", pencil, Y, X) ;
   exit(-1) ;
}
if ( PENCIL_IS_COMPLEX(pencil) && !DENSEMTX_IS_COMPLEX(X) ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n pencil is complex, X is not\n", pencil, Y, X) ;
   exit(-1) ;
}
DenseMtx_dimensions(Y, &nrowY, &ncolY) ;
DenseMtx_dimensions(X, &nrowX, &ncolX) ;
if ( nrowY != nrowX || ncolY != ncolX ) {
   fprintf(stderr, "\n fatal error in Pencil_mmm(%p,%p,%p)"
           "\n nrowY %d, ncolY %d, nrowX %d, ncolX %d\n", 
           pencil, Y, X, nrowY, ncolY, nrowX, ncolX) ;
   exit(-1) ;
}
nrow = nrowY ;
nrhs = ncolY ;
if ( pencil->inpmtxA == NULL ) {
/*
   -----------------
   A is the identity
   -----------------
*/
   if ( PENCIL_IS_REAL(pencil) ) {
      double   *x, *y ;
      int      irow, jrhs ;

      x = DenseMtx_entries(X) ;
      y = DenseMtx_entries(Y) ;
      for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            y[irow] += x[irow] ;
         }
         x += nrow ; y += nrow ;
      }
   } else if ( PENCIL_IS_COMPLEX(pencil) ) {
      double   *x, *y ;
      int      irow, jrhs ;

      for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
         for ( irow = 0 ; irow < nrow ; irow++ ) {
            y[2*irow]   += x[2*irow]   ;
            y[2*irow+1] += x[2*irow+1] ;
         }
         x += 2*nrow ; y += 2*nrow ;
      }
   }
} else {
   double   alpha[2] ;
/*
   ----------------------------------------
   A is not the identity, multiply x with A
   ----------------------------------------
*/
   alpha[0] = 1.0 ; alpha[1] = 0.0 ;
   if ( PENCIL_IS_SYMMETRIC(pencil) ) {
      InpMtx_sym_mmm(pencil->inpmtxA, Y, alpha, X) ;
   } else if ( PENCIL_IS_HERMITIAN(pencil) ) {
      InpMtx_herm_mmm(pencil->inpmtxA, Y, alpha, X) ;
   } else if ( PENCIL_IS_NONSYMMETRIC(pencil) ) {
      InpMtx_nonsym_mmm(pencil->inpmtxA, Y, alpha, X) ;
   }
}
if ( pencil->sigma[0] != 0.0 || pencil->sigma[1] != 0.0 ) {
   if ( pencil->inpmtxB != NULL ) {
/*
      -----------------------------------------
      B is not the identity, add sigma*B*x to y
      -----------------------------------------
*/
      if ( PENCIL_IS_SYMMETRIC(pencil) ) {
         InpMtx_sym_mmm(pencil->inpmtxB, Y, pencil->sigma, X) ;
      } else if ( PENCIL_IS_HERMITIAN(pencil) ) {
         InpMtx_herm_mmm(pencil->inpmtxB, Y, pencil->sigma, X) ;
      } else if ( PENCIL_IS_NONSYMMETRIC(pencil) ) {
         InpMtx_nonsym_mmm(pencil->inpmtxB, Y, pencil->sigma, X) ;
      }
   } else {
/*
      -----------------------------------
      B is the identity, add sigma*x to y
      -----------------------------------
*/
      if ( PENCIL_IS_REAL(pencil) ) {
         double   sigmareal ;
         double   *x, *y ;
         int      irow, jrhs ;

         x = DenseMtx_entries(X) ; y = DenseMtx_entries(Y) ;
         sigmareal = pencil->sigma[0] ;
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nrow ; irow++ ) {
               y[irow] += sigmareal*x[irow] ;
            }
            x += nrow ; y += nrow ;
         }
      } else if ( PENCIL_IS_COMPLEX(pencil) ) {
         double   sigmaimag, sigmareal, ximag, xreal ;
         double   *x, *y ;
         int      irow, jrhs ;

         x = DenseMtx_entries(X) ; y = DenseMtx_entries(Y) ;
         sigmareal = pencil->sigma[0] ; sigmaimag = pencil->sigma[1] ;
         for ( jrhs = 0 ; jrhs < nrhs ; jrhs++ ) {
            for ( irow = 0 ; irow < nrow ; irow++ ) {
               xreal = x[2*irow] ; ximag = x[2*irow+1] ;
               y[2*irow]   += sigmareal*xreal - sigmaimag*ximag ;
               y[2*irow+1] += sigmareal*ximag + sigmaimag*xreal ;
            }
            x += 2*nrow ; y += 2*nrow ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
