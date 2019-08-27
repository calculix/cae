/*  loadEntries.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   load entries from sigma*A

   chv     -- pointer to the Chv object that holds the front
   pencil  -- pointer to a Pencil that holds the matrix entries 
   msglvl  -- message level
   msgFile -- message file

   created  -- 97jul18, cca
   ------------------------------------------------------------
*/
void
FrontMtx_loadEntries (
   Chv      *chv,
   Pencil   *pencil,
   int      msglvl,
   FILE     *msgFile
) {
InpMtx   *inpmtxA, *inpmtxB ;
double   one[2] = {1.0,0.0} ;
double   *sigma ;
double   *chvent ;
int      chvsize, ichv, ncol, nD, nL, nU ;
int      *chvind, *colind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_loadEntries(%p,%p,%d,%p)"
           "\n bad input\n", chv, pencil, msglvl, msgFile) ;
   exit(-1) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, 
           "\n\n # inside loadEntries for chv %d" 
           ", sigma = %12.4e + i*%12.4e",
           chv->id, pencil->sigma[0], pencil->sigma[1]) ;
   fflush(msgFile) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
Chv_columnIndices(chv, &ncol, &colind) ;
/*
   ----------------------------------------
   load the original entries, A + sigma * B
   ----------------------------------------
*/
inpmtxA = pencil->inpmtxA ;
sigma   = pencil->sigma   ;
inpmtxB = pencil->inpmtxB ;
if ( inpmtxA != NULL ) {
   int   ii ;
/*
   -------------------
   load entries from A
   -------------------
*/
   for ( ii = 0 ; ii < nD ; ii++ ) {
      ichv = colind[ii] ;
      if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) { 
         InpMtx_realVector(inpmtxA, ichv, &chvsize, &chvind, &chvent) ;
      } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) { 
         InpMtx_complexVector(inpmtxA, 
                              ichv, &chvsize, &chvind, &chvent) ;
      }
      if ( chvsize > 0 ) {
         if ( msglvl > 3 ) {
            int ierr ;
            fprintf(msgFile, "\n inpmtxA chevron %d : chvsize = %d", 
                    ichv, chvsize) ;
            fprintf(msgFile, "\n chvind") ;
            IVfp80(msgFile, chvsize, chvind, 80, &ierr) ;
            fprintf(msgFile, "\n chvent") ;
            if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) { 
               DVfprintf(msgFile, chvsize, chvent) ;
            } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) { 
               DVfprintf(msgFile, 2*chvsize, chvent) ;
            }
            fflush(msgFile) ;
         }
         Chv_addChevron(chv, one, ichv, chvsize, chvind, chvent) ;
      }
   }
} else {
   double   *entries ;
   int      ii, off, stride ;
/*
   -----------------
   load the identity
   -----------------
*/
   entries = Chv_entries(chv) ;
   if ( CHV_IS_REAL(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         stride = nD + chv->nU ;
         off    = 0 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[off] += 1.0 ;
            off += stride ;
            stride-- ;
         }
      } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
         stride = 2*nD + chv->nL + chv->nU - 2 ;
         off    = nD + chv->nL - 1 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[off] += 1.0 ;
            off += stride ;
            stride -= 2 ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         stride = nD + chv->nU ;
         off    = 0 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[2*off] += 1.0 ;
            off += stride ;
            stride-- ;
         }
      } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
         stride = 2*nD + chv->nL + chv->nU - 2 ;
         off    = nD + chv->nL - 1 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[2*off] += 1.0 ;
            off += stride ;
            stride -= 2 ;
         }
      }
   }
}
if ( inpmtxB != NULL ) {
   int   ii ;
/*
   -------------------------
   load entries from sigma*B
   -------------------------
*/
   for ( ii = 0 ; ii < nD ; ii++ ) {
      ichv = colind[ii] ;
      if ( INPMTX_IS_REAL_ENTRIES(inpmtxB) ) { 
         InpMtx_realVector(inpmtxB, ichv, &chvsize, &chvind, &chvent) ;
      } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) { 
         InpMtx_complexVector(inpmtxB, 
                              ichv, &chvsize, &chvind, &chvent) ;
      }
      if ( chvsize > 0 ) {
         if ( msglvl > 3 ) {
            int ierr ;
            fprintf(msgFile, "\n inpmtxB chevron %d : chvsize = %d", 
                    ichv, chvsize) ;
            fprintf(msgFile, "\n chvind") ;
            IVfp80(msgFile, chvsize, chvind, 80, &ierr) ;
            fprintf(msgFile, "\n chvent") ;
            if ( INPMTX_IS_REAL_ENTRIES(inpmtxA) ) { 
               DVfprintf(msgFile, chvsize, chvent) ;
            } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtxA) ) { 
               DVfprintf(msgFile, 2*chvsize, chvent) ;
            }
         }
         Chv_addChevron(chv, sigma, ichv, chvsize, chvind, chvent) ;
      }
   }
} else {
   double   *entries ;
   int      ii, off, stride ;
/*
   --------------------------------------
   load a scalar multiple of the identity
   --------------------------------------
*/
   entries = Chv_entries(chv) ;
   if ( CHV_IS_REAL(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) ) {
         stride = nD + chv->nU ;
         off    = 0 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[off] += sigma[0] ;
            off += stride ;
            stride-- ;
         }
      } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
         stride = 2*nD + chv->nL + chv->nU - 2 ;
         off    = nD + chv->nL - 1 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[off] += sigma[0] ;
            off += stride ;
            stride -= 2 ;
         }
      }
   } else if ( CHV_IS_COMPLEX(chv) ) {
      if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
         if ( CHV_IS_HERMITIAN(chv) && sigma[1] != 0.0 ) {
            fprintf(stderr, 
                    "\n fatal error in FrontMtx_loadEntries()"
                    "\n chevron is hermitian" 
                    "\n sigma = %12.4e + %12.4e*i\n",
                    sigma[0], sigma[1]) ;
            exit(-1) ;
         }
         stride = nD + chv->nU ;
         off    = 0 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[2*off]   += sigma[0] ;
            entries[2*off+1] += sigma[1] ;
            off += stride ;
            stride-- ;
         }
      } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
         stride = 2*nD + chv->nL + chv->nU - 2 ;
         off    = nD + chv->nL - 1 ;
         for ( ii = 0 ; ii < nD ; ii++ ) {
            entries[2*off]   += sigma[0] ;
            entries[2*off+1] += sigma[1] ;
            off += stride ;
            stride -= 2 ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
