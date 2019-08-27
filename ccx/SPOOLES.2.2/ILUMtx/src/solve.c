/*  solve.c  */

#include "../ILUMtx.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- to solve the linear system
      (L + I)D(I + U) X = B, (U^T + I)D(I + U) X = B or
      (U^H + I)D(I + U) X = B. X and B are single vectors.

   workDV is a temporary vector. 
   note, if workDV is different than B, then B is unchanged on return.
   one can have X, B and workDV all point to the same object.
 
   if pops is not NULL, then on return *pops has been incremented
   by the number of operations performed in the solve.

   return values ---
      1  -- normal return
     -1  -- mtx is NULL
     -2  -- neqns <= 0
     -3  -- bad type for mtx
     -4  -- bad symmetryflag for mtx
     -5  -- storage mode of L is invalid
     -6  -- storage mode of U is invalid
     -7  -- sizesL is NULL
     -8  -- sizesU is NULL
     -9  -- p_indL is NULL
     -10 -- p_indU is NULL
     -11 -- entD is NULL
     -12 -- p_entL is NULL
     -13 -- p_entU is NULL
     -14 -- X is NULL
     -15 -- size of X is incorrect
     -16 -- entries of X are NULL
     -17 -- B is NULL
     -18 -- size of B is incorrect
     -19 -- entries of B are NULL
     -20 -- workDV is NULL
     -21 -- size of workDV != neqns
     -22 -- entries of workDV are NULL
     -23 -- msglvl > 0 and msgFile is NULL

   created -- 98oct03, cca
   -------------------------------------------------------------------
*/
int
ILUMtx_solveVector (
   ILUMtx   *mtx,
   DV       *X,
   DV       *B,
   DV       *workDV,
   double   *pops,
   int      msglvl,
   FILE     *msgFile
) {
double   ops ;
double   *bent, *entD, *work, *xent ;
double   **p_entL, **p_entU ;
int      bsize, ieqn, neqns, wsize, xsize ;
int      *sizesL, *sizesU ;
int      **p_indL, **p_indU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), mtx = NULL\n") ;
   return(-1) ;
}
if ( (neqns = mtx->neqns) <= 0 ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n neqns = %d\n", neqns) ;
   return(-2) ;
}
if ( !(ILUMTX_IS_REAL(mtx) || ILUMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n type = %d\n", mtx->type) ;
   return(-3) ;
}
if ( !(ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx)
       || ILUMTX_IS_NONSYMMETRIC(mtx)) ) {
   fprintf(stderr, "\n error in ILUMsolveVector()"
           "\n mtx symmetry = %d\n", mtx->symmetryflag) ;
   return(-4) ;
}
if ( !(ILUMTX_IS_L_BY_ROWS(mtx) || ILUMTX_IS_L_BY_COLUMNS(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n LstorageMode = %d\n", mtx->LstorageMode) ;
   return(-5) ;
}
if ( !(ILUMTX_IS_U_BY_ROWS(mtx) || ILUMTX_IS_U_BY_COLUMNS(mtx)) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n UstorageMode = %d\n", mtx->UstorageMode) ;
   return(-6) ;
}
sizesU = mtx->sizesU ;
entD   = mtx->entD   ;
p_indU = mtx->p_indU ;
p_entU = mtx->p_entU ;
if ( ILUMTX_IS_SYMMETRIC(mtx) || ILUMTX_IS_HERMITIAN(mtx) ) {
   sizesL = mtx->sizesU ;
   p_indL = mtx->p_indU ;
   p_entL = mtx->p_entU ;
} else {
   sizesL = mtx->sizesL ;
   p_indL = mtx->p_indL ;
   p_entL = mtx->p_entL ;
}
if ( sizesL == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), sizesL = NULL\n") ;
   return(-7) ;
}
if ( sizesU == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), sizesU = NULL\n") ;
   return(-8) ;
}
if ( p_indL == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), p_indL = NULL\n") ;
   return(-9) ;
}
if ( p_indU == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), p_indU = NULL\n") ;
   return(-10) ;
}
if ( entD == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), entD = NULL\n") ;
   return(-11) ;
}
if ( p_entL == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), p_entL = NULL\n") ;
   return(-12) ;
}
if ( p_entU == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), p_entU = NULL\n") ;
   return(-13) ;
}
if ( X == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), X = NULL\n") ;
   return(-14) ;
}
DV_sizeAndEntries(X, &xsize, &xent) ;
if (  (ILUMTX_IS_REAL(mtx) && xsize != neqns)
   || (ILUMTX_IS_COMPLEX(mtx) && xsize != 2*neqns) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n neqns = %d, size of X = %d\n", neqns, xsize) ;
   return(-15) ;
}
if ( xent == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n entries of X are NULL\n") ;
   return(-16) ;
}
if ( B == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), B = NULL\n") ;
   return(-17) ;
}
DV_sizeAndEntries(B, &bsize, &bent) ;
if (  (ILUMTX_IS_REAL(mtx) && bsize != neqns)
   || (ILUMTX_IS_COMPLEX(mtx) && bsize != 2*neqns) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n neqns = %d, size of B = %d\n", neqns, bsize) ;
   return(-18) ;
}
if ( bent == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n entries of B are NULL\n") ;
   return(-19) ;
}
if ( workDV == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector(), workDV = NULL\n") ;
   return(-20) ;
}
DV_sizeAndEntries(workDV, &wsize, &work) ;
if (  (ILUMTX_IS_REAL(mtx) && wsize != neqns)
   || (ILUMTX_IS_COMPLEX(mtx) && wsize != 2*neqns) ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n neqns = %d, size of workDV = %d\n", neqns, wsize) ;
   return(-21) ;
}
if ( work == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n entries of workDV are NULL\n") ;
   return(-22) ;
}
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n error in ILUM_solveVector()"
           "\n msglvl = %d, msgFile is NULL\n", msglvl) ;
   return(-23) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------
   copy B vector into work[]
   -------------------------
*/
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n bent") ;
   DVfprintf(msgFile, neqns, bent) ;
}
if ( bent != work ) {
   if ( ILUMTX_IS_REAL(mtx) ) {
      DVcopy(neqns, work, bent) ;
   } else {
      DVcopy(2*neqns, work, bent) ;
   }
}
/*
   -------------
   solve L Y = B 
   -------------
*/
ops = 0.0 ;
if ( ILUMTX_IS_REAL(mtx) ) {
   if ( mtx->LstorageMode == SPOOLES_BY_COLUMNS ) {
      int      ii, Lsize, *ind ;
      double   yi, *ent ;
   
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (Lsize = sizesL[ieqn]) > 0 ) {
            yi  = work[ieqn]   ;
            ind = p_indL[ieqn] ;
            ent = p_entL[ieqn] ;
            for ( ii = 0 ; ii < Lsize ; ii++ ) {
               work[ind[ii]] -= yi * ent[ii] ;
            }
            ops += 2*Lsize ;
         }
      }
   } else {
      int      ii, Lsize, *ind ;
      double   sum, *ent ;
   
      for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
         if ( (Lsize = sizesL[ieqn]) > 0 ) {
            ind = p_indL[ieqn] ;
            ent = p_entL[ieqn] ;
            for ( ii = 0, sum = 0.0 ; ii < Lsize ; ii++ ) {
               sum += work[ind[ii]] * ent[ii] ;
            }
            work[ieqn] -= sum ;
            ops += 2*Lsize ;
         }
      }
   }
} else {
   if ( mtx->LstorageMode == SPOOLES_BY_COLUMNS ) {
      int      ii, jj, kk, Lsize, *ind ;
      double   LI, LR, yI, yR, *ent ;
   
      if ( ILUMTX_IS_HERMITIAN(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (Lsize = sizesL[ieqn]) > 0 ) {
               yR  = work[2*ieqn]   ;
               yI  = work[2*ieqn+1] ;
               ind = p_indL[ieqn]   ;
               ent = p_entL[ieqn]   ;
               for ( ii = kk = 0 ; ii < Lsize ; ii++, kk += 2 ) {
                  LR = ent[kk] ; LI = -ent[kk+1] ;
                  jj   = 2*ind[ii] ; 
                  work[jj]   -= yR*LR - yI*LI ;
                  work[jj+1] -= yR*LI + yI*LR ;
               }
               ops += 8*Lsize ;
            }
         }
      } else {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (Lsize = sizesL[ieqn]) > 0 ) {
               yR  = work[2*ieqn]   ;
               yI  = work[2*ieqn+1] ;
               ind = p_indL[ieqn]   ;
               ent = p_entL[ieqn]   ;
               for ( ii = kk = 0 ; ii < Lsize ; ii++, kk += 2 ) {
                  LR = ent[kk] ; LI = ent[kk+1] ;
                  jj   = 2*ind[ii] ; 
                  work[jj]   -= yR*LR - yI*LI ;
                  work[jj+1] -= yR*LI + yI*LR ;
               }
               ops += 8*Lsize ;
            }
         }
      }
   } else {
      int      ii, jj, kk, Lsize, *ind ;
      double   LI, LR, sumI, sumR, yI, yR, *ent ;
   
      if ( ILUMTX_IS_HERMITIAN(mtx) ) {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (Lsize = sizesL[ieqn]) > 0 ) {
               ind  = p_indL[ieqn] ;
               ent  = p_entL[ieqn] ;
               sumI = sumR = 0.0 ;
               for ( ii = kk = 0 ; ii < Lsize ; ii++, kk += 2 ) {
                  jj = 2*ind[ii] ; yR = work[jj] ; yI = work[jj+1] ;
                  LR = ent[kk] ; LI = -ent[kk+1] ;
                  sumR += yR*LR - yI*LI ;
                  sumI += yR*LI + yI*LR ;
               }
               work[2*ieqn]   -= sumR ;
               work[2*ieqn+1] -= sumI ;
               ops += 8*Lsize ;
            }
         }
      } else {
         for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
            if ( (Lsize = sizesL[ieqn]) > 0 ) {
               ind  = p_indL[ieqn] ;
               ent  = p_entL[ieqn] ;
               sumI = sumR = 0.0 ;
               for ( ii = kk = 0 ; ii < Lsize ; ii++, kk += 2 ) {
                  jj = 2*ind[ii] ; yR = work[jj] ; yI = work[jj+1] ;
                  LR = ent[kk] ; LI = ent[kk+1] ;
                  sumR += yR*LR - yI*LI ;
                  sumI += yR*LI + yI*LR ;
               }
               work[2*ieqn]   -= sumR ;
               work[2*ieqn+1] -= sumI ;
               ops += 8*Lsize ;
            }
         }
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %% after forward solve") ;
   fprintf(msgFile, "\n Y = zeros(%d,1) ;", neqns) ;
   if ( ILUMTX_IS_REAL(mtx) ) {
      DV_writeForMatlab(workDV, "Y", msgFile) ;
   } else {
      workDV->size = neqns ;
      ZV_writeForMatlab((ZV *) workDV, "Y", msgFile) ;
      workDV->size = 2*neqns ;
   }
   fflush(msgFile) ;
}
/*
   -------------
   solve D Z = Y 
   -------------
*/
if ( ILUMTX_IS_REAL(mtx) ) {
   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      work[ieqn] = work[ieqn] / entD[ieqn] ;
   }
   ops += 2*neqns ;
} else {
   double   dI, dR, yI, yR ;
   int      rc ;

   for ( ieqn = 0 ; ieqn < neqns ; ieqn++ ) {
      rc = Zrecip(entD[2*ieqn], entD[2*ieqn+1], &dR, &dI) ;
      if ( rc == 1 ) {
         yR = work[2*ieqn] ; yI = work[2*ieqn+1] ;
         work[2*ieqn]   = yR*dR - yI*dI ;
         work[2*ieqn+1] = yR*dI + yI*dR ;
      }
   }
   ops += 11*neqns ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %% after diagonal solve") ;
   fprintf(msgFile, "\n Z = zeros(%d,1) ;", neqns) ;
   if ( ILUMTX_IS_REAL(mtx) ) {
      DV_writeForMatlab(workDV, "Z", msgFile) ;
   } else {
      workDV->size = neqns ;
      ZV_writeForMatlab((ZV *) workDV, "Z", msgFile) ;
      workDV->size = 2*neqns ;
   }
   fflush(msgFile) ;
}
/*
   -------------
   solve U X = Z 
   -------------
*/
if ( ILUMTX_IS_REAL(mtx) ) {
   if ( mtx->UstorageMode == SPOOLES_BY_COLUMNS ) {
      int      ii, Usize, *ind ;
      double   xi, *ent ;
   
      for ( ieqn = neqns - 1 ; ieqn >= 0 ; ieqn-- ) {
         if ( (Usize = sizesU[ieqn]) > 0 ) {
            xi  = work[ieqn]   ;
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0 ; ii < Usize ; ii++ ) {
               work[ind[ii]] -= xi * ent[ii] ;
            }
            ops += 2*Usize ;
         }
      }
   } else {
      int      ii, Usize, *ind ;
      double   sum, *ent ;
   
      for ( ieqn = neqns - 1 ; ieqn >= 0 ; ieqn-- ) {
         if ( (Usize = sizesU[ieqn]) > 0 ) {
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = 0, sum = 0.0 ; ii < Usize ; ii++ ) {
               sum += work[ind[ii]] * ent[ii] ;
            }
            work[ieqn] -= sum ;
            ops += 2*Usize ;
         }
      }
   }
} else {
   if ( mtx->UstorageMode == SPOOLES_BY_COLUMNS ) {
      int      ii, jj, kk, Usize, *ind ;
      double   UI, UR, xI, xR, *ent ;
   
      for ( ieqn = neqns - 1 ; ieqn >= 0 ; ieqn-- ) {
         if ( (Usize = sizesU[ieqn]) > 0 ) {
            xR  = work[2*ieqn] ; xI = work[2*ieqn+1] ;
            ind = p_indU[ieqn] ;
            ent = p_entU[ieqn] ;
            for ( ii = kk = 0 ; ii < Usize ; ii++, kk += 2 ) {
               jj = 2*ind[ii] ;
               UR = ent[kk] ; UI = ent[kk+1] ;
               work[jj]   -= xR*UR - xI*UI ;
               work[jj+1] -= xR*UI + xI*UR ;
            }
            ops += 8*Usize ;
         }
      }
   } else {
      int      ii, jj, kk, Usize, *ind ;
      double   sumI, sumR, UI, UR, xI, xR, *ent ;
   
      for ( ieqn = neqns - 1 ; ieqn >= 0 ; ieqn-- ) {
         if ( (Usize = sizesU[ieqn]) > 0 ) {
            ind  = p_indU[ieqn] ;
            ent  = p_entU[ieqn] ;
            sumR = sumI = 0.0 ;
            for ( ii = kk = 0 ; ii < Usize ; ii++, kk += 2 ) {
               jj = 2*ind[ii] ; xR = work[jj] ; xI = work[jj+1] ;
               UR = ent[kk] ; UI = ent[kk+1] ;
               sumR += xR*UR - xI*UI ;
               sumI += xR*UI + xI*UR ;
            }
            work[2*ieqn]   -= sumR ;
            work[2*ieqn+1] -= sumI ;
            ops += 8*Usize ;
         }
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %% after backward solve") ;
   fprintf(msgFile, "\n W = zeros(%d,1) ;", neqns) ;
   if ( ILUMTX_IS_REAL(mtx) ) {
      DV_writeForMatlab(workDV, "W", msgFile) ;
   } else {
      workDV->size = neqns ;
      ZV_writeForMatlab((ZV *) workDV, "W", msgFile) ;
      workDV->size = 2*neqns ;
   }
   fflush(msgFile) ;
}
/*
   -----------------------
   copy work vector into X
   -----------------------
*/
if ( work != xent ) {
   if ( ILUMTX_IS_REAL(mtx) ) {
      DVcopy(neqns, xent, work) ;
   } else {
      DVcopy(2*neqns, xent, work) ;
   }
}
/*
   -----------------------------------------
   increment the operation count if not NULL
   -----------------------------------------
*/
if ( pops != NULL ) {
   *pops += ops ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
