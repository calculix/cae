/*  update.c  */

#include "../Chv.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose --  perform the hermitian factor update 
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= U_{I, \bnd{I} \cap J}^H D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= U_{I, \bnd{I} \cap J}^H D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}

   created -- 98apr17, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateH (
   Chv      *chvT,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV       *tempDV
) {
int   firstInT, firstInU, jcolT, jcolU, lastInT, lastInU, ncolT, ncolU ;
int   *colindT, *colindU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvT == NULL || mtxD == NULL || mtxU == NULL || tempDV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_updateH(%p,%p,%p,%p)"
           "\n bad input\n", chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
if ( ! CHV_IS_COMPLEX(chvT) ) {
   fprintf(stderr, "\n fatal error in Chv_updateH(%p,%p,%p,%p)"
           "\n bad input, chvT is not complex\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtxD) ) {
   fprintf(stderr, "\n fatal error in Chv_updateH(%p,%p,%p,%p)"
           "\n bad input, mtxD is not complex\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtxU) ) {
   fprintf(stderr, "\n fatal error in Chv_updateH(%p,%p,%p,%p)"
           "\n bad input, mtxU is not complex\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
Chv_columnIndices(chvT, &ncolT, &colindT) ;
SubMtx_columnIndices(mtxU, &ncolU, &colindU) ;
/*
   -----------------------------
   locate first column of U in T
   -----------------------------
*/
firstInT = colindT[0] ;
lastInT  = colindT[chvT->nD-1] ;
for ( jcolU = 0 ; jcolU < ncolU ; jcolU++ ) {
   if ( firstInT <= colindU[jcolU] && colindU[jcolU] <= lastInT ) {
      break ;
   }
}
if ( (firstInU = jcolU) == ncolU ) {
   return ;
}
/*
   ----------------------------
   locate last column of U in T
   ----------------------------
*/
lastInU = firstInU ;
for (    ; jcolU < ncolU ; jcolU++ ) {
   if ( colindU[jcolU] <= lastInT ) {
      lastInU = jcolU ;
   } else {
      break ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% firstInU = %d, lastInU = %d", 
        firstInU, lastInU) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------------------
   overwrite supported column indices with indices local to T
   ----------------------------------------------------------
*/
for ( jcolU = firstInU, jcolT = 0 ; jcolU < ncolU ; jcolU++ ) {
   while ( colindU[jcolU] != colindT[jcolT] ) {
      jcolT++ ;
   }
   colindU[jcolU] = jcolT ;
}
if ( SUBMTX_IS_DENSE_COLUMNS(mtxU) ) {
   double   isum, rsum ;
   double   sums[18] ;
   double   *base0, *base1, *base2, *colU0, *colU1, *colU2, *entU,
            *rowUT0, *rowUT1, *rowUT2, *temp0, *temp1, *temp2 ;
   int      ichv0, ichv1, ichv2, ii, inc1, inc2, irowUT, 
            kloc0, kloc1, kloc2, nrowU ;

   SubMtx_denseInfo(mtxU, &nrowU, &ncolU, &inc1, &inc2, &entU) ;
   DV_setSize(tempDV, 6*nrowU) ;
   temp0 = DV_entries(tempDV) ;
   temp1 = temp0 + 2*nrowU ;
   temp2 = temp1 + 2*nrowU ;
/*
   --------------------------------------------
   loop over the rows of U^H in groups of three
   --------------------------------------------
*/
   rowUT0 = entU + 2*firstInU*nrowU ;
   for ( irowUT = firstInU ; irowUT <= lastInU - 2 ; irowUT += 3 ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% working on rows %d, %d and %d",
              colindU[irowUT], colindU[irowUT+1], colindU[irowUT+2]) ;
      fflush(stdout) ;
#endif
      rowUT1 = rowUT0 + 2*nrowU ;
      rowUT2 = rowUT1 + 2*nrowU ;
/*
      -----------------------------------------------------
      get base locations for the chevron's diagonal entries
      -----------------------------------------------------
*/
      ichv0 = colindU[irowUT] ;
      base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
      ichv1 = colindU[irowUT+1] ;
      base1 = Chv_diagLocation(chvT, ichv1) - 2*ichv1 ;
      ichv2 = colindU[irowUT+2] ;
      base2 = Chv_diagLocation(chvT, ichv2) - 2*ichv2 ;
/*
      ------------------------------------
              [ temp0 ]   [ rowUT0 ]^H
      compute [ temp1 ] = [ rowUT1 ]   * D
              [ temp2 ]   [ rowUT2 ]
      ------------------------------------
*/
      DVzero(6*nrowU, temp0) ;
      SubMtx_scale3vec(mtxD, temp0, temp1, temp2,
                       rowUT0, rowUT1, rowUT2) ;
      for ( ii = 0 ; ii < nrowU ; ii++ ) {
         temp0[2*ii+1] = -temp0[2*ii+1] ;
         temp1[2*ii+1] = -temp1[2*ii+1] ;
         temp2[2*ii+1] = -temp2[2*ii+1] ;
      }
/*
      --------------------------------------------------
      update the 3x3 upper triangle for these three rows
      --------------------------------------------------
*/
      colU0 = rowUT0 ;
      colU1 = rowUT1 ;
      colU2 = rowUT2 ;
      ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base0[2*ichv0+1]) ;
}
*/
/*
      base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
*/
      base0[2*ichv0] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e",
           isum, base0[2*ichv0+1]) ;
}
*/
      ZVdotU(nrowU, temp0, colU1, &rsum, &isum) ;
      base0[2*ichv1] -= rsum ; base0[2*ichv1+1] -= isum ;
      ZVdotU(nrowU, temp0, colU2, &rsum, &isum) ;
      base0[2*ichv2] -= rsum ; base0[2*ichv2+1] -= isum ;
      ZVdotU(nrowU, temp1, colU1, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base1[2*ichv1+1]) ;
}
*/
/*
      base1[2*ichv1] -= rsum ; base1[2*ichv1+1] -= isum ;
*/
      base1[2*ichv1] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e",
           isum, base1[2*ichv1+1]) ;
}
*/
      ZVdotU(nrowU, temp1, colU2, &rsum, &isum) ;
      base1[2*ichv2] -= rsum ; base1[2*ichv2+1] -= isum ;
      ZVdotU(nrowU, temp2, colU2, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base2[2*ichv2+1]) ;
}
*/
/*
      base2[2*ichv2] -= rsum ; base2[2*ichv2+1] -= isum ;
*/
      base2[2*ichv2] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e",
           isum, base2[2*ichv2+1]) ;
}
*/
      colU0 = colU2 + 2*nrowU ;
/*
      --------------------------------------
      update the remainder of the three rows
      --------------------------------------
*/
      for ( jcolU = irowUT + 3 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
         colU1 = colU0 + 2*nrowU ;
         colU2 = colU1 + 2*nrowU ;
         ZVdotU33(nrowU, temp0, temp1, temp2, 
                  colU0, colU1, colU2, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         kloc2 = 2*colindU[jcolU+2] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
         base1[kloc0] -= sums[ 6] ; base1[kloc0+1] -= sums[ 7] ;
         base1[kloc1] -= sums[ 8] ; base1[kloc1+1] -= sums[ 9] ;
         base1[kloc2] -= sums[10] ; base1[kloc2+1] -= sums[11] ;
         base2[kloc0] -= sums[12] ; base2[kloc0+1] -= sums[13] ;
         base2[kloc1] -= sums[14] ; base2[kloc1+1] -= sums[15] ;
         base2[kloc2] -= sums[16] ; base2[kloc2+1] -= sums[17] ;
         colU0 = colU2 + 2*nrowU ;
      }
      if ( jcolU == ncolU - 2 ) {
         colU1 = colU0 + 2*nrowU ;
         ZVdotU32(nrowU, temp0, temp1, temp2, colU0, colU1, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         base1[kloc0] -= sums[ 4] ; base1[kloc0+1] -= sums[ 5] ;
         base1[kloc1] -= sums[ 6] ; base1[kloc1+1] -= sums[ 7] ;
         base2[kloc0] -= sums[ 8] ; base2[kloc0+1] -= sums[ 9] ;
         base2[kloc1] -= sums[10] ; base2[kloc1+1] -= sums[11] ;
      } else if ( jcolU == ncolU - 1 ) {
         ZVdotU31(nrowU, temp0, temp1, temp2, colU0, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base1[kloc0] -= sums[ 2] ; base1[kloc0+1] -= sums[ 3] ;
         base2[kloc0] -= sums[ 4] ; base2[kloc0+1] -= sums[ 5] ;
      }
      rowUT0 = rowUT2 + 2*nrowU ;
   }
   if ( irowUT == lastInU - 1 ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% working on rows %d and %d",
              colindU[irowUT], colindU[irowUT+1]) ;
      fflush(stdout) ;
#endif
      rowUT1 = rowUT0 + 2*nrowU ;
/*
      -----------------------------------------------------
      get base locations for the chevron's diagonal entries
      -----------------------------------------------------
*/
      ichv0 = colindU[irowUT] ;
      base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
      ichv1 = colindU[irowUT+1] ;
      base1 = Chv_diagLocation(chvT, ichv1) - 2*ichv1 ;
/*
      ------------------------------------
              [ temp0 ]   [ rowUT0 ]^H
      compute [ temp1 ] = [ rowUT1 ]   * D
      ------------------------------------
*/
      DVzero(4*nrowU, temp0) ;
      SubMtx_scale2vec(mtxD, temp0, temp1, rowUT0, rowUT1) ;
      for ( ii = 0 ; ii < nrowU ; ii++ ) {
         temp0[2*ii+1] = -temp0[2*ii+1] ;
         temp1[2*ii+1] = -temp1[2*ii+1] ;
      }
/*
      ------------------------------------------------
      update the 2x2 upper triangle for these two rows
      ------------------------------------------------
*/
      colU0 = rowUT0 ;
      colU1 = rowUT1 ;
      ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base0[2*ichv0+1]) ;
}
*/
/*
      base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
*/
      base0[2*ichv0] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e", 
           isum, base0[2*ichv0+1]) ;
}
*/
      ZVdotU(nrowU, temp0, colU1, &rsum, &isum) ;
      base0[2*ichv1] -= rsum ; base0[2*ichv1+1] -= isum ;
      ZVdotU(nrowU, temp1, colU1, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base1[2*ichv1+1]) ;
}
*/
/*
      base1[2*ichv1] -= rsum ; base1[2*ichv1+1] -= isum ;
*/
      base1[2*ichv1] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e", 
           isum, base1[2*ichv1+1]) ;
}
*/
      colU0 = colU1 + 2*nrowU ;
/*
      ------------------------------------
      update the remainder of the two rows
      ------------------------------------
*/
      for ( jcolU = irowUT + 2 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
         colU1 = colU0 + 2*nrowU ;
         colU2 = colU1 + 2*nrowU ;
         ZVdotU23(nrowU, temp0, temp1, colU0, colU1, colU2, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         kloc2 = 2*colindU[jcolU+2] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
         base1[kloc0] -= sums[ 6] ; base1[kloc0+1] -= sums[ 7] ;
         base1[kloc1] -= sums[ 8] ; base1[kloc1+1] -= sums[ 9] ;
         base1[kloc2] -= sums[10] ; base1[kloc2+1] -= sums[11] ;
         colU0 = colU2 + 2*nrowU ;
      }
      if ( jcolU == ncolU - 2 ) {
         colU1 = colU0 + 2*nrowU ;
         ZVdotU22(nrowU, temp0, temp1, colU0, colU1, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         base1[kloc0] -= sums[ 4] ; base1[kloc0+1] -= sums[ 5] ;
         base1[kloc1] -= sums[ 6] ; base1[kloc1+1] -= sums[ 7] ;
      } else if ( jcolU == ncolU - 1 ) {
         ZVdotU21(nrowU, temp0, temp1, colU0, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base1[kloc0] -= sums[ 2] ; base1[kloc0+1] -= sums[ 3] ;
      }
   } else if ( irowUT == lastInU ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n %% working on row %d", colindU[irowUT]) ;
      fflush(stdout) ;
#endif
/*
      -----------------------------------------------------
      get base locations for the chevron's diagonal entries
      -----------------------------------------------------
*/
      ichv0 = colindU[irowUT] ;
      base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
/*
      ------------------------------------
      compute [ temp0 ] = [ rowUT0 ]^H * D
      ------------------------------------
*/
      DVzero(2*nrowU, temp0) ;
      SubMtx_scale1vec(mtxD, temp0, rowUT0) ;
      for ( ii = 0 ; ii < nrowU ; ii++ ) {
         temp0[2*ii+1] = -temp0[2*ii+1] ;
      }
/*
      ------------------------------------------
      update the 1x1 upper triangle for this row
      ------------------------------------------
*/
      colU0 = rowUT0 ;
      ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, "\n (%d,%d) : imag(diag) = %12.5e",
           chvT->id, mtxD->rowid, base0[2*ichv0+1]) ;
}
*/
/*
      base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
*/
      base0[2*ichv0] -= rsum ; 
/*
if ( chvT->id != -1 ) {
   fprintf(stdout, ", isum = %12.5e, imag(diag) = %12.5e", 
           isum, base0[2*ichv0+1]) ;
}
*/
      colU0 = colU0 + 2*nrowU ;
/*
      -------------------------------
      update the remainder of the row
      -------------------------------
*/
      for ( jcolU = irowUT + 1 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on columns %d, %d and %d", 
                 jcolU, jcolU+1, jcolU+2) ;
         fflush(stdout) ;
#endif
         colU1 = colU0 + 2*nrowU ;
         colU2 = colU1 + 2*nrowU ;
         ZVdotU13(nrowU, temp0, colU0, colU1, colU2, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         kloc2 = 2*colindU[jcolU+2] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
         colU0 = colU2 + 2*nrowU ;
      }
      if ( jcolU == ncolU - 2 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on columns %d and %d", 
                 jcolU, jcolU+1) ;
         fflush(stdout) ;
#endif
         colU1 = colU0 + 2*nrowU ;
         ZVdotU12(nrowU, temp0, colU0, colU1, sums) ;
         kloc0 = 2*colindU[jcolU] ;
         kloc1 = 2*colindU[jcolU+1] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
      } else if ( jcolU == ncolU - 1 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on column %d", jcolU) ;
         fflush(stdout) ;
#endif
         ZVdotU11(nrowU, temp0, colU0, sums) ;
/*
fprintf(stdout, "\n SUMS[0] = %12.4e, SUMS[1] = %12.4e",
        sums[0], sums[1]) ;
*/
         kloc0 = 2*colindU[jcolU] ;
         base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
/*
fprintf(stdout, "\n base0[%d] = %12.4e, base0[%d] = %12.4e",
        kloc0, base0[kloc0], kloc0+1, base0[kloc0+1]) ;
*/
      }
   }
} else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxU) ) {
   double   isum, rsum ;
   double   *base0, *colU0, *entU, *rowUT0, *temp0, *temp1 ;
   int      ichv0, ii, iloc, irowUT, kloc0, nentU, nrowU, offset, 
            rloc, sizeU, sizeUT ;
   int      *indU, *indU0, *indUT0, *sizes ;

   SubMtx_sparseColumnsInfo(mtxU, &ncolU, &nentU, &sizes, &indU, &entU) ;
   nrowU = mtxU->nrow ;
   DV_setSize(tempDV, 4*nrowU) ;
   temp0 = DV_entries(tempDV) ;
   temp1 = temp0 + 2*nrowU ;
/*
   -------------------------------------------
   get the offset into the indices and entries
   -------------------------------------------
*/
   for ( jcolU = offset = 0 ; jcolU < firstInU ; jcolU++ ) {
      offset += sizes[jcolU] ;
   }
/*
   ------------------------------------
   loop over the supporting rows of U^H
   ------------------------------------
*/
   rowUT0 = entU + 2*offset ;
   indUT0 = indU + offset ;
   for ( irowUT = firstInU ; irowUT <= lastInU ; irowUT++ ) {
      if ( (sizeUT = sizes[irowUT]) > 0 ) {
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
/*
         ------------------------------------
         compute [ temp0 ] = [ rowUT0 ]^H * D
         ------------------------------------
*/
         DVzero(4*nrowU, temp0) ;
         for ( ii = 0 ; ii < sizeUT ; ii++ ) {
            rloc = 2*indUT0[ii] ; iloc = rloc + 1 ;
            temp1[rloc] = rowUT0[2*ii] ;
            temp1[iloc] = rowUT0[2*ii+1] ;
         }
         SubMtx_scale1vec(mtxD, temp0, temp1) ;
         for ( ii = 0 ; ii < nrowU ; ii++ ) {
            temp0[2*ii+1] = -temp0[2*ii+1] ;
         }
/*
         -------------------------------
         loop over the following columns
         -------------------------------
*/
         colU0 = rowUT0 ;
         indU0 = indUT0 ;
         for ( jcolU = irowUT ; jcolU < ncolU ; jcolU++ ) {
            if ( (sizeU = sizes[jcolU]) > 0 ) {
               ZVdotiU(sizeU, temp0, indU0, colU0, &rsum, &isum) ;
               kloc0 = 2*colindU[jcolU] ;
               base0[kloc0]   -= rsum ; base0[kloc0+1] -= isum ;
               colU0 += 2*sizeU ;
               indU0 += sizeU ;
            }
         }
         rowUT0 += 2*sizeUT ;
         indUT0 += sizeUT ;
      }
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_updateH(%p,%p,%p,%p)"
           "\n mtxU must have dense or sparse columns\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
/*
   ---------------------------------------------------
   overwrite the local indices with the global indices
   ---------------------------------------------------
*/
for ( jcolU = firstInU ; jcolU < ncolU ; jcolU++ ) {
   colindU[jcolU] = colindT[colindU[jcolU]] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose --  perform the symmetric factor update 
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= U_{I, \bnd{I} \cap J}^T D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= U_{I, \bnd{I} \cap J}^T D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}

   created -- 98apr17, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateS (
   Chv      *chvT,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV       *tempDV
) {
int   firstInT, firstInU, jcolT, jcolU, lastInT, lastInU, ncolT, ncolU ;
int   *colindT, *colindU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvT == NULL || mtxD == NULL || mtxU == NULL || tempDV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_updateS(%p,%p,%p,%p)"
           "\n bad input\n", chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chvT) ) {
   if ( ! SUBMTX_IS_REAL(mtxD) || ! SUBMTX_IS_REAL(mtxU) ) {
      fprintf(stderr, "\n fatal error in Chv_updateT(%p,%p,%p,%p)"
              "\n chvT is real, but mtxD and/or mtxU are not\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( ! SUBMTX_IS_COMPLEX(mtxD) || ! SUBMTX_IS_COMPLEX(mtxU) ) {
      fprintf(stderr, "\n fatal error in Chv_updateT(%p,%p,%p,%p)"
              "\n chvT is complex, but mtxD and/or mtxU are not\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_updateT(%p,%p,%p,%p)"
           "\n bad input, chvT is not real or complex\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
Chv_columnIndices(chvT, &ncolT, &colindT) ;
SubMtx_columnIndices(mtxU, &ncolU, &colindU) ;
#if MYDEBUG > 0
fprintf(stdout, "\n %% Chv column indices") ;
IVfprintf(stdout, ncolT, colindT) ;
fprintf(stdout, "\n %% submtx column indices") ;
IVfprintf(stdout, ncolU, colindU) ;
fflush(stdout) ;
#endif
/*
   -----------------------------
   locate first column of U in T
   -----------------------------
*/
firstInT = colindT[0] ;
lastInT  = colindT[chvT->nD-1] ;
for ( jcolU = 0 ; jcolU < ncolU ; jcolU++ ) {
   if ( firstInT <= colindU[jcolU] && colindU[jcolU] <= lastInT ) {
      break ;
   }
}
if ( (firstInU = jcolU) == ncolU ) {
   return ;
}
/*
   ----------------------------
   locate last column of U in T
   ----------------------------
*/
lastInU = firstInU ;
for (    ; jcolU < ncolU ; jcolU++ ) {
   if ( colindU[jcolU] <= lastInT ) {
      lastInU = jcolU ;
   } else {
      break ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% firstInU = %d, lastInU = %d", 
        firstInU, lastInU) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------------------
   overwrite supported column indices with indices local to T
   ----------------------------------------------------------
*/
for ( jcolU = firstInU, jcolT = 0 ; jcolU < ncolU ; jcolU++ ) {
   while ( colindU[jcolU] != colindT[jcolT] ) {
      jcolT++ ;
   }
   colindU[jcolU] = jcolT ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% submtx column indices") ;
IVfprintf(stdout, ncolU, colindU) ;
fflush(stdout) ;
#endif
if ( CHV_IS_REAL(chvT) ) {
   if ( SUBMTX_IS_DENSE_COLUMNS(mtxU) ) {
      double   sums[9] ;
      double   *base0, *base1, *base2, *colU0, *colU1, *colU2, *entU,
               *rowUT0, *rowUT1, *rowUT2, *temp0, *temp1, *temp2 ;
      int      ichv0, ichv1, ichv2, inc1, inc2, irowUT, 
               kloc0, kloc1, kloc2, nrowU ;

      SubMtx_denseInfo(mtxU, &nrowU, &ncolU, &inc1, &inc2, &entU) ;
      DV_setSize(tempDV, 3*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + nrowU ;
      temp2 = temp1 + nrowU ;
/*
      --------------------------------------------
      loop over the rows of U^T in groups of three
      --------------------------------------------
*/
      rowUT0 = entU + firstInU*nrowU ;
      for ( irowUT = firstInU ; irowUT <= lastInU - 2 ; irowUT += 3 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d, %d and %d",
              colindU[irowUT], colindU[irowUT+1], colindU[irowUT+2]) ;
         fflush(stdout) ;
#endif
         rowUT1 = rowUT0 + nrowU ;
         rowUT2 = rowUT1 + nrowU ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - ichv0 ;
         ichv1 = colindU[irowUT+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) - ichv1 ;
         ichv2 = colindU[irowUT+2] ;
         base2 = Chv_diagLocation(chvT, ichv2) - ichv2 ;
/*
         ------------------------------------
                 [ temp0 ]   [ rowUT0 ]^T
         compute [ temp1 ] = [ rowUT1 ]   * D
                 [ temp2 ]   [ rowUT2 ]
         ------------------------------------
*/
         DVzero(3*nrowU, temp0) ;
         SubMtx_scale3vec(mtxD, temp0, temp1, temp2,
                          rowUT0, rowUT1, rowUT2) ;
/*
         --------------------------------------------------
         update the 3x3 upper triangle for these three rows
         --------------------------------------------------
*/
         colU0 = rowUT0 ;
         colU1 = rowUT1 ;
         colU2 = rowUT2 ;
         base0[ichv0] -= DVdot(nrowU, temp0, colU0) ;
         base0[ichv1] -= DVdot(nrowU, temp0, colU1) ;
         base0[ichv2] -= DVdot(nrowU, temp0, colU2) ;
         base1[ichv1] -= DVdot(nrowU, temp1, colU1) ;
         base1[ichv2] -= DVdot(nrowU, temp1, colU2) ;
         base2[ichv2] -= DVdot(nrowU, temp2, colU2) ;
         colU0 = colU2 + nrowU ;
/*
         --------------------------------------
         update the remainder of the three rows
         --------------------------------------
*/
         for ( jcolU = irowUT + 3 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
            DVdot33(nrowU, temp0, temp1, temp2, 
                     colU0, colU1, colU2, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            kloc2 = colindU[jcolU+2] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
            base0[kloc2] -= sums[2] ; 
            base1[kloc0] -= sums[3] ; 
            base1[kloc1] -= sums[4] ; 
            base1[kloc2] -= sums[5] ; 
            base2[kloc0] -= sums[6] ; 
            base2[kloc1] -= sums[7] ; 
            base2[kloc2] -= sums[8] ; 
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            colU1 = colU0 + nrowU ;
            DVdot32(nrowU, temp0, temp1, temp2, colU0, colU1, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
            base1[kloc0] -= sums[2] ; 
            base1[kloc1] -= sums[3] ; 
            base2[kloc0] -= sums[4] ; 
            base2[kloc1] -= sums[5] ; 
         } else if ( jcolU == ncolU - 1 ) {
            DVdot31(nrowU, temp0, temp1, temp2, colU0, sums) ;
            kloc0 = colindU[jcolU] ;
            base0[kloc0] -= sums[0] ; 
            base1[kloc0] -= sums[1] ; 
            base2[kloc0] -= sums[2] ; 
         }
         rowUT0 = rowUT2 + nrowU ;
      }
      if ( irowUT == lastInU - 1 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d and %d",
                 colindU[irowUT], colindU[irowUT+1]) ;
         fflush(stdout) ;
#endif
         rowUT1 = rowUT0 + nrowU ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - ichv0 ;
         ichv1 = colindU[irowUT+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) - ichv1 ;
/*
         ------------------------------------
                 [ temp0 ]   [ rowUT0 ]^T
         compute [ temp1 ] = [ rowUT1 ]   * D
         ------------------------------------
*/
         DVzero(2*nrowU, temp0) ;
         SubMtx_scale2vec(mtxD, temp0, temp1, rowUT0, rowUT1) ;
/*
         ------------------------------------------------
         update the 2x2 upper triangle for these two rows
         ------------------------------------------------
*/
         colU0 = rowUT0 ;
         colU1 = rowUT1 ;
         base0[ichv0] -= DVdot(nrowU, temp0, colU0) ;
         base0[ichv1] -= DVdot(nrowU, temp0, colU1) ;
         base1[ichv1] -= DVdot(nrowU, temp1, colU1) ;
         colU0 = colU1 + nrowU ;
/*
         ------------------------------------
         update the remainder of the two rows
         ------------------------------------
*/
         for ( jcolU = irowUT + 2 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
            DVdot23(nrowU, temp0, temp1, colU0, colU1, colU2, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            kloc2 = colindU[jcolU+2] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
            base0[kloc2] -= sums[2] ; 
            base1[kloc0] -= sums[3] ; 
            base1[kloc1] -= sums[4] ; 
            base1[kloc2] -= sums[5] ; 
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            colU1 = colU0 + nrowU ;
            DVdot22(nrowU, temp0, temp1, colU0, colU1, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
            base1[kloc0] -= sums[2] ; 
            base1[kloc1] -= sums[3] ; 
         } else if ( jcolU == ncolU - 1 ) {
            DVdot21(nrowU, temp0, temp1, colU0, sums) ;
            kloc0 = colindU[jcolU] ;
            base0[kloc0] -= sums[0] ; 
            base1[kloc0] -= sums[1] ; 
         }
      } else if ( irowUT == lastInU ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on row %d", colindU[irowUT]) ;
         fflush(stdout) ;
#endif
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - ichv0 ;
/*
         ------------------------------------
         compute [ temp0 ] = [ rowUT0 ]^T * D
         ------------------------------------
*/
         DVzero(nrowU, temp0) ;
         SubMtx_scale1vec(mtxD, temp0, rowUT0) ;
/*
         ------------------------------------------
         update the 1x1 upper triangle for this row
         ------------------------------------------
*/
         colU0 = rowUT0 ;
         base0[ichv0] -= DVdot(nrowU, temp0, colU0) ;
         colU0 = colU0 + nrowU ;
/*
         -------------------------------
         update the remainder of the row
         -------------------------------
*/
         for ( jcolU = irowUT + 1 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on columns %d, %d and %d", 
                    jcolU, jcolU+1, jcolU+2) ;
            fflush(stdout) ;
#endif
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
            DVdot13(nrowU, temp0, colU0, colU1, colU2, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            kloc2 = colindU[jcolU+2] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
            base0[kloc2] -= sums[2] ; 
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on columns %d and %d", 
                    jcolU, jcolU+1) ;
            fflush(stdout) ;
#endif
            colU1 = colU0 + nrowU ;
            DVdot12(nrowU, temp0, colU0, colU1, sums) ;
            kloc0 = colindU[jcolU] ;
            kloc1 = colindU[jcolU+1] ;
            base0[kloc0] -= sums[0] ; 
            base0[kloc1] -= sums[1] ; 
         } else if ( jcolU == ncolU - 1 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on column %d", jcolU) ;
            fflush(stdout) ;
#endif
            DVdot11(nrowU, temp0, colU0, sums) ;
            kloc0 = colindU[jcolU] ;
            base0[kloc0] -= sums[0] ; 
         }
      }
   } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxU) ) {
      double   sum ;
      double   *base0, *colU0, *entU, *rowUT0, *temp0, *temp1 ;
      int      ichv0, ii, irowUT, kloc0, loc, nentU, nrowU, offset, 
               sizeU, sizeUT ;
      int      *indU, *indU0, *indUT0, *sizes ;
   
      SubMtx_sparseColumnsInfo(mtxU, 
                               &ncolU, &nentU, &sizes, &indU, &entU) ;
      nrowU = mtxU->nrow ;
      DV_setSize(tempDV, 2*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + nrowU ;
/*
      -------------------------------------------
      get the offset into the indices and entries
      -------------------------------------------
*/
      for ( jcolU = offset = 0 ; jcolU < firstInU ; jcolU++ ) {
         offset += sizes[jcolU] ;
      }
/*
      ------------------------------------
      loop over the supporting rows of U^T
      ------------------------------------
*/
      rowUT0 = entU + offset ;
      indUT0 = indU + offset ;
      for ( irowUT = firstInU ; irowUT <= lastInU ; irowUT++ ) {
         if ( (sizeUT = sizes[irowUT]) > 0 ) {
/*
            -----------------------------------------------------
            get base locations for the chevron's diagonal entries
            -----------------------------------------------------
*/
            ichv0 = colindU[irowUT] ;
            base0 = Chv_diagLocation(chvT, ichv0) - ichv0 ;
/*
            ------------------------------------
            compute [ temp0 ] = [ rowUT0 ]^T * D
            ------------------------------------
*/
            DVzero(2*nrowU, temp0) ;
            for ( ii = 0 ; ii < sizeUT ; ii++ ) {
               loc = indUT0[ii] ; 
               temp1[loc] = rowUT0[ii] ;
            }
            SubMtx_scale1vec(mtxD, temp0, temp1) ;
/*
            -------------------------------
            loop over the following columns
            -------------------------------
*/
            colU0 = rowUT0 ;
            indU0 = indUT0 ;
            for ( jcolU = irowUT ; jcolU < ncolU ; jcolU++ ) {
               if ( (sizeU = sizes[jcolU]) > 0 ) {
                  sum = DVdoti(sizeU, temp0, indU0, colU0) ;
                  kloc0 = colindU[jcolU] ;
                  base0[kloc0] -= sum ; 
                  colU0 += sizeU ;
                  indU0 += sizeU ;
               }
            }
            rowUT0 += sizeUT ;
            indUT0 += sizeUT ;
         }
      }
   } else {
      fprintf(stderr, "\n fatal error in Chv_updateT(%p,%p,%p,%p)"
              "\n bad input, mtxU must have dense or sparse columns\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( SUBMTX_IS_DENSE_COLUMNS(mtxU) ) {
      double   isum, rsum ;
      double   sums[18] ;
      double   *base0, *base1, *base2, *colU0, *colU1, *colU2, *entU,
               *rowUT0, *rowUT1, *rowUT2, *temp0, *temp1, *temp2 ;
      int      ichv0, ichv1, ichv2, inc1, inc2, irowUT, 
               kloc0, kloc1, kloc2, nrowU ;

      SubMtx_denseInfo(mtxU, &nrowU, &ncolU, &inc1, &inc2, &entU) ;
      DV_setSize(tempDV, 6*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + 2*nrowU ;
      temp2 = temp1 + 2*nrowU ;
/*
      --------------------------------------------
      loop over the rows of U^T in groups of three
      --------------------------------------------
*/
      rowUT0 = entU + 2*firstInU*nrowU ;
      for ( irowUT = firstInU ; irowUT <= lastInU - 2 ; irowUT += 3 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d, %d and %d",
              colindU[irowUT], colindU[irowUT+1], colindU[irowUT+2]) ;
         fflush(stdout) ;
#endif
         rowUT1 = rowUT0 + 2*nrowU ;
         rowUT2 = rowUT1 + 2*nrowU ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
         ichv1 = colindU[irowUT+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) - 2*ichv1 ;
         ichv2 = colindU[irowUT+2] ;
         base2 = Chv_diagLocation(chvT, ichv2) - 2*ichv2 ;
/*
         ------------------------------------
                 [ temp0 ]   [ rowUT0 ]^T
         compute [ temp1 ] = [ rowUT1 ]   * D
                 [ temp2 ]   [ rowUT2 ]
         ------------------------------------
*/
         DVzero(6*nrowU, temp0) ;
         SubMtx_scale3vec(mtxD, temp0, temp1, temp2,
                          rowUT0, rowUT1, rowUT2) ;
/*
         --------------------------------------------------
         update the 3x3 upper triangle for these three rows
         --------------------------------------------------
*/
         colU0 = rowUT0 ;
         colU1 = rowUT1 ;
         colU2 = rowUT2 ;
         ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
         base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
         ZVdotU(nrowU, temp0, colU1, &rsum, &isum) ;
         base0[2*ichv1] -= rsum ; base0[2*ichv1+1] -= isum ;
         ZVdotU(nrowU, temp0, colU2, &rsum, &isum) ;
         base0[2*ichv2] -= rsum ; base0[2*ichv2+1] -= isum ;
         ZVdotU(nrowU, temp1, colU1, &rsum, &isum) ;
         base1[2*ichv1] -= rsum ; base1[2*ichv1+1] -= isum ;
         ZVdotU(nrowU, temp1, colU2, &rsum, &isum) ;
         base1[2*ichv2] -= rsum ; base1[2*ichv2+1] -= isum ;
         ZVdotU(nrowU, temp2, colU2, &rsum, &isum) ;
         base2[2*ichv2] -= rsum ; base2[2*ichv2+1] -= isum ;
         colU0 = colU2 + 2*nrowU ;
/*
         --------------------------------------
         update the remainder of the three rows
         --------------------------------------
*/
         for ( jcolU = irowUT + 3 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
            ZVdotU33(nrowU, temp0, temp1, temp2, 
                     colU0, colU1, colU2, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            kloc2 = 2*colindU[jcolU+2] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
            base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
            base1[kloc0] -= sums[ 6] ; base1[kloc0+1] -= sums[ 7] ;
            base1[kloc1] -= sums[ 8] ; base1[kloc1+1] -= sums[ 9] ;
            base1[kloc2] -= sums[10] ; base1[kloc2+1] -= sums[11] ;
            base2[kloc0] -= sums[12] ; base2[kloc0+1] -= sums[13] ;
            base2[kloc1] -= sums[14] ; base2[kloc1+1] -= sums[15] ;
            base2[kloc2] -= sums[16] ; base2[kloc2+1] -= sums[17] ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            colU1 = colU0 + 2*nrowU ;
            ZVdotU32(nrowU, temp0, temp1, temp2, colU0, colU1, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
            base1[kloc0] -= sums[ 4] ; base1[kloc0+1] -= sums[ 5] ;
            base1[kloc1] -= sums[ 6] ; base1[kloc1+1] -= sums[ 7] ;
            base2[kloc0] -= sums[ 8] ; base2[kloc0+1] -= sums[ 9] ;
            base2[kloc1] -= sums[10] ; base2[kloc1+1] -= sums[11] ;
         } else if ( jcolU == ncolU - 1 ) {
            ZVdotU31(nrowU, temp0, temp1, temp2, colU0, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base1[kloc0] -= sums[ 2] ; base1[kloc0+1] -= sums[ 3] ;
            base2[kloc0] -= sums[ 4] ; base2[kloc0+1] -= sums[ 5] ;
         }
         rowUT0 = rowUT2 + 2*nrowU ;
      }
      if ( irowUT == lastInU - 1 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d and %d",
                 colindU[irowUT], colindU[irowUT+1]) ;
         fflush(stdout) ;
#endif
         rowUT1 = rowUT0 + 2*nrowU ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
         ichv1 = colindU[irowUT+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) - 2*ichv1 ;
/*
         ------------------------------------
                 [ temp0 ]   [ rowUT0 ]^T
         compute [ temp1 ] = [ rowUT1 ]   * D
         ------------------------------------
*/
         DVzero(4*nrowU, temp0) ;
         SubMtx_scale2vec(mtxD, temp0, temp1, rowUT0, rowUT1) ;
/*
         ------------------------------------------------
         update the 2x2 upper triangle for these two rows
         ------------------------------------------------
*/
         colU0 = rowUT0 ;
         colU1 = rowUT1 ;
         ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
         base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
         ZVdotU(nrowU, temp0, colU1, &rsum, &isum) ;
         base0[2*ichv1] -= rsum ; base0[2*ichv1+1] -= isum ;
         ZVdotU(nrowU, temp1, colU1, &rsum, &isum) ;
         base1[2*ichv1] -= rsum ; base1[2*ichv1+1] -= isum ;
         colU0 = colU1 + 2*nrowU ;
/*
         ------------------------------------
         update the remainder of the two rows
         ------------------------------------
*/
         for ( jcolU = irowUT + 2 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
            ZVdotU23(nrowU, temp0, temp1, colU0, colU1, colU2, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            kloc2 = 2*colindU[jcolU+2] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
            base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
            base1[kloc0] -= sums[ 6] ; base1[kloc0+1] -= sums[ 7] ;
            base1[kloc1] -= sums[ 8] ; base1[kloc1+1] -= sums[ 9] ;
            base1[kloc2] -= sums[10] ; base1[kloc2+1] -= sums[11] ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            colU1 = colU0 + 2*nrowU ;
            ZVdotU22(nrowU, temp0, temp1, colU0, colU1, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
            base1[kloc0] -= sums[ 4] ; base1[kloc0+1] -= sums[ 5] ;
            base1[kloc1] -= sums[ 6] ; base1[kloc1+1] -= sums[ 7] ;
         } else if ( jcolU == ncolU - 1 ) {
            ZVdotU21(nrowU, temp0, temp1, colU0, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base1[kloc0] -= sums[ 2] ; base1[kloc0+1] -= sums[ 3] ;
         }
      } else if ( irowUT == lastInU ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on row %d", colindU[irowUT]) ;
         fflush(stdout) ;
#endif
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowUT] ;
         base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
/*
         ------------------------------------
         compute [ temp0 ] = [ rowUT0 ]^T * D
         ------------------------------------
*/
         DVzero(2*nrowU, temp0) ;
         SubMtx_scale1vec(mtxD, temp0, rowUT0) ;
/*
         ------------------------------------------
         update the 1x1 upper triangle for this row
         ------------------------------------------
*/
         colU0 = rowUT0 ;
         ZVdotU(nrowU, temp0, colU0, &rsum, &isum) ;
         base0[2*ichv0] -= rsum ; base0[2*ichv0+1] -= isum ;
         colU0 = colU0 + 2*nrowU ;
/*
         -------------------------------
         update the remainder of the row
         -------------------------------
*/
         for ( jcolU = irowUT + 1 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on columns %d, %d and %d", 
                    jcolU, jcolU+1, jcolU+2) ;
            fflush(stdout) ;
#endif
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
            ZVdotU13(nrowU, temp0, colU0, colU1, colU2, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            kloc2 = 2*colindU[jcolU+2] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
            base0[kloc2] -= sums[ 4] ; base0[kloc2+1] -= sums[ 5] ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on columns %d and %d", 
                    jcolU, jcolU+1) ;
            fflush(stdout) ;
#endif
            colU1 = colU0 + 2*nrowU ;
            ZVdotU12(nrowU, temp0, colU0, colU1, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            kloc1 = 2*colindU[jcolU+1] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
            base0[kloc1] -= sums[ 2] ; base0[kloc1+1] -= sums[ 3] ;
         } else if ( jcolU == ncolU - 1 ) {
#if MYDEBUG > 0
            fprintf(stdout, "\n %% working on column %d", jcolU) ;
            fflush(stdout) ;
#endif
            ZVdotU11(nrowU, temp0, colU0, sums) ;
            kloc0 = 2*colindU[jcolU] ;
            base0[kloc0] -= sums[ 0] ; base0[kloc0+1] -= sums[ 1] ;
         }
      }
   } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxU) ) {
      double   isum, rsum ;
      double   *base0, *colU0, *entU, *rowUT0, *temp0, *temp1 ;
      int      ichv0, ii, iloc, irowUT, kloc0, nentU, nrowU, offset, 
               rloc, sizeU, sizeUT ;
      int      *indU, *indU0, *indUT0, *sizes ;
   
      SubMtx_sparseColumnsInfo(mtxU, 
                               &ncolU, &nentU, &sizes, &indU, &entU) ;
      nrowU = mtxU->nrow ;
      DV_setSize(tempDV, 4*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + 2*nrowU ;
/*
      -------------------------------------------
      get the offset into the indices and entries
      -------------------------------------------
*/
      for ( jcolU = offset = 0 ; jcolU < firstInU ; jcolU++ ) {
         offset += sizes[jcolU] ;
      }
/*
      ------------------------------------
      loop over the supporting rows of U^T
      ------------------------------------
*/
      rowUT0 = entU + 2*offset ;
      indUT0 = indU + offset ;
      for ( irowUT = firstInU ; irowUT <= lastInU ; irowUT++ ) {
         if ( (sizeUT = sizes[irowUT]) > 0 ) {
/*
            -----------------------------------------------------
            get base locations for the chevron's diagonal entries
            -----------------------------------------------------
*/
            ichv0 = colindU[irowUT] ;
            base0 = Chv_diagLocation(chvT, ichv0) - 2*ichv0 ;
/*
            ------------------------------------
            compute [ temp0 ] = [ rowUT0 ]^T * D
            ------------------------------------
*/
            DVzero(4*nrowU, temp0) ;
            for ( ii = 0 ; ii < sizeUT ; ii++ ) {
               rloc = 2*indUT0[ii] ; iloc = rloc + 1 ;
               temp1[rloc] = rowUT0[2*ii] ;
               temp1[iloc] = rowUT0[2*ii+1] ;
            }
            SubMtx_scale1vec(mtxD, temp0, temp1) ;
/*
            -------------------------------
            loop over the following columns
            -------------------------------
*/
            colU0 = rowUT0 ;
            indU0 = indUT0 ;
            for ( jcolU = irowUT ; jcolU < ncolU ; jcolU++ ) {
               if ( (sizeU = sizes[jcolU]) > 0 ) {
                  ZVdotiU(sizeU, temp0, indU0, colU0, &rsum, &isum) ;
                  kloc0 = 2*colindU[jcolU] ;
                  base0[kloc0]   -= rsum ; base0[kloc0+1] -= isum ;
                  colU0 += 2*sizeU ;
                  indU0 += sizeU ;
               }
            }
            rowUT0 += 2*sizeUT ;
            indUT0 += sizeUT ;
         }
      }
   } else {
      fprintf(stderr, "\n fatal error in Chv_updateT(%p,%p,%p,%p)"
              "\n bad input, mtxU must have dense or sparse columns\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
}
/*
   ---------------------------------------------------
   overwrite the local indices with the global indices
   ---------------------------------------------------
*/
for ( jcolU = firstInU ; jcolU < ncolU ; jcolU++ ) {
   colindU[jcolU] = colindT[colindU[jcolU]] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   purpose --  perform the nonsymmetric factor update 
     T_{\bnd{I} \cap J, \bnd{I} \cap J}
          -= L_{\bnd{I} \cap J, I} D_{I, I} U_{I, \bnd{I} \cap J}
   and
     T_{\bnd{I} \cap J, \bnd{I} \cap \bnd{J}}
         -= L_{\bnd{I} \cap J, I} D_{I, I} U_{I, \bnd{I} \cap \bnd{J}}
   and
     T_{\bnd{I} \cap \bnd{J}, \bnd{I} \cap J}
         -= L_{\bnd{I} \cap \bnd{J}, I} D_{I, I} U_{I, \bnd{I} \cap J}

   created -- 98feb27, cca
   ---------------------------------------------------------------------
*/
void
Chv_updateN (
   Chv      *chvT,
   SubMtx   *mtxL,
   SubMtx   *mtxD,
   SubMtx   *mtxU,
   DV       *tempDV
) {
int   firstInT, firstInU, jcolT, jcolU, lastInT, lastInU, ncolT, ncolU ;
int   *colindT, *colindU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chvT == NULL || mtxD == NULL || mtxU == NULL || tempDV == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p,%p)"
           "\n bad input\n", chvT, mtxL, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chvT) ) {
   if ( ! SUBMTX_IS_REAL(mtxD) || ! SUBMTX_IS_REAL(mtxL) 
     || ! SUBMTX_IS_REAL(mtxU) ) {
      fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p)"
              "\n chvT is real, but mtxD, mtxL and/or mtxU are not\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( ! SUBMTX_IS_COMPLEX(mtxD) || ! SUBMTX_IS_COMPLEX(mtxL) 
     || ! SUBMTX_IS_COMPLEX(mtxU) ) {
      fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p)"
             "\n chvT is complex, but mtxD, mtxL and/or mtxU are not\n",
             chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p)"
           "\n bad input, chvT is not real or complex\n", 
           chvT, mtxD, mtxU, tempDV) ;
   exit(-1) ;
}
Chv_columnIndices(chvT, &ncolT, &colindT) ;
SubMtx_columnIndices(mtxU, &ncolU, &colindU) ;
/*
   -----------------------------
   locate first column of U in T
   -----------------------------
*/
firstInT = colindT[0] ;
lastInT  = colindT[chvT->nD-1] ;
for ( jcolU = 0 ; jcolU < ncolU ; jcolU++ ) {
   if ( firstInT <= colindU[jcolU] && colindU[jcolU] <= lastInT ) {
      break ;
   }
}
if ( (firstInU = jcolU) == ncolU ) {
   return ;
}
/*
   ----------------------------
   locate last column of U in T
   ----------------------------
*/
lastInU = firstInU ;
for (    ; jcolU < ncolU ; jcolU++ ) {
   if ( colindU[jcolU] <= lastInT ) {
      lastInU = jcolU ;
   } else {
      break ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %% firstInU = %d, lastInU = %d", 
        firstInU, lastInU) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------------------
   overwrite supported column indices with indices local to T
   ----------------------------------------------------------
*/
for ( jcolU = firstInU, jcolT = 0 ; jcolU < ncolU ; jcolU++ ) {
   while ( colindU[jcolU] != colindT[jcolT] ) {
      jcolT++ ;
   }
   colindU[jcolU] = jcolT ;
}
if ( CHV_IS_REAL(chvT) ) {
   if ( SUBMTX_IS_DENSE_COLUMNS(mtxU) && SUBMTX_IS_DENSE_ROWS(mtxL) ) {
      double   sums[9] ;
      double   *base0, *base1, *base2, *colU0, *colU1, *colU2, *entL,
               *entU, *Ltemp0, *Ltemp1, *Ltemp2, *rowL0, *rowL1, *rowL2,
               *Utemp0, *Utemp1, *Utemp2 ;
      int      ichv0, ichv1, ichv2, inc1, inc2, irowL, jcolU, 
               loc, loc0, loc1, loc2, ncolL, nrowL, nrowU, offset ;

      SubMtx_denseInfo(mtxL, &nrowL, &ncolL, &inc1, &inc2, &entL) ;
      SubMtx_denseInfo(mtxU, &nrowU, &ncolU, &inc1, &inc2, &entU) ;
      DV_setSize(tempDV, 6*nrowU) ;
      Ltemp0 = DV_entries(tempDV) ;
      Ltemp1 = Ltemp0 + nrowU ;
      Ltemp2 = Ltemp1 + nrowU ;
      Utemp0 = Ltemp2 + nrowU ;
      Utemp1 = Utemp0 + nrowU ;
      Utemp2 = Utemp1 + nrowU ;
/*
      -----------------------------------------------------
      loop over the rows of L in groups of three
      updating the diagonal and upper blocks of the chevron
      -----------------------------------------------------
*/
      offset = firstInU*nrowU ;
      for ( irowL = firstInU ; irowL <= lastInU - 2 ; irowL += 3 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d, %d and %d",
                 colindU[irowL], colindU[irowL+1], colindU[irowL+2]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + offset ;
         rowL1 = rowL0 + nrowU  ;
         rowL2 = rowL1 + nrowU  ;
         colU0 = entU  + offset ;
         colU1 = colU0 + nrowU  ;
         colU2 = colU1 + nrowU  ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
         ichv1 = colindU[irowL+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) ;
         ichv2 = colindU[irowL+2] ;
         base2 = Chv_diagLocation(chvT, ichv2) ;
/*
         -----------------------------------
                 [ Ltemp0 ]   [ rowL0 ]
         compute [ Ltemp1 ] = [ rowL1 ] * D
                 [ Ltemp2 ]   [ rowL2 ]
         -----------------------------------
*/
         DVzero(3*nrowU, Ltemp0) ;
         SubMtx_scale3vec(mtxD, Ltemp0, Ltemp1, Ltemp2,
                          rowL0, rowL1, rowL2) ;
/*
         -----------------------------------
                 [ Utemp0 ]       [ colU0 ]
         compute [ Utemp1 ] = D * [ colU0 ]
                 [ Utemp2 ]       [ colU0 ]
         -----------------------------------
*/
         DVzero(3*nrowU, Utemp0) ;
         SubMtx_scale3vec(mtxD, Utemp0, Utemp1, Utemp2,
                          colU0, colU1, colU2) ;
/*
         --------------------------------------------------------------
         update the 3x3 diagonal block for these three rows and columns
         --------------------------------------------------------------
*/
         DVdot33(nrowU, Ltemp0, Ltemp1, Ltemp2, 
                 colU0, colU1, colU2, sums) ;
/*
         -------------------------------------
         inject the nine sums into the chevron
         -------------------------------------
*/
         base0[0]    -= sums[0] ; 
         loc = ichv1 - ichv0 ;
         base0[loc]  -= sums[1] ; 
         base0[-loc] -= sums[3] ; 
         loc = ichv2 - ichv0 ;
         base0[loc]  -= sums[2] ; 
         base0[-loc] -= sums[6] ; 
         base1[0]    -= sums[4] ; 
         loc = ichv2 - ichv1 ;
         base1[loc]  -= sums[5] ; 
         base1[-loc] -= sums[7] ; 
         base2[0]    -= sums[8] ; 
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL2 + nrowU ;
         colU0 = colU2 + nrowU ;
         for ( jcolU = irowL + 3 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + nrowU ;
            rowL2 = rowL1 + nrowU ;
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
            loc2 = colindU[jcolU+2] ;
/*
            ---------------------
            upper 3x3 block first
            ---------------------
*/
            DVdot33(nrowU, Ltemp0, Ltemp1, Ltemp2, 
                    colU0, colU1, colU2, sums) ;
/*
            -------------------------------------
            inject the nine sums into the chevron
            -------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base0[loc2] -= sums[2] ; 
            base1[loc0] -= sums[3] ; 
            base1[loc1] -= sums[4] ; 
            base1[loc2] -= sums[5] ; 
            base2[loc0] -= sums[6] ; 
            base2[loc1] -= sums[7] ; 
            base2[loc2] -= sums[8] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
/*
            ----------------------
            lower 3x3 block second
            ----------------------
*/
            DVdot33(nrowU, rowL0, rowL1, rowL2, 
                    Utemp0, Utemp1, Utemp2, sums) ;
/*
            -------------------------------------
            inject the nine sums into the chevron
            -------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base2[-loc0] -= sums[2] ; 
            base0[-loc1] -= sums[3] ; 
            base1[-loc1] -= sums[4] ; 
            base2[-loc1] -= sums[5] ; 
            base0[-loc2] -= sums[6] ; 
            base1[-loc2] -= sums[7] ; 
            base2[-loc2] -= sums[8] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + nrowU ;
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + nrowU ;
            colU1 = colU0 + nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
/*
            ---------------------
            upper 3x2 block first
            ---------------------
*/
            DVdot32(nrowU, Ltemp0, Ltemp1, Ltemp2, colU0, colU1, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base1[loc0] -= sums[2] ; 
            base1[loc1] -= sums[3] ; 
            base2[loc0] -= sums[4] ; 
            base2[loc1] -= sums[5] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
/*
            ----------------------
            lower 2x3 block second
            ----------------------
*/
            DVdot23(nrowU, rowL0, rowL1, Utemp0, Utemp1, Utemp2, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base2[-loc0] -= sums[2] ; 
            base0[-loc1] -= sums[3] ; 
            base1[-loc1] -= sums[4] ; 
            base2[-loc1] -= sums[5] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = colindU[jcolU] ;
/*
            ---------------------
            upper 3x1 block first
            ---------------------
*/
            DVdot31(nrowU, Ltemp0, Ltemp1, Ltemp2, colU0, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
            base0[loc0] -= sums[0] ; 
            base1[loc0] -= sums[1] ; 
            base2[loc0] -= sums[2] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
/*
            ----------------------
            lower 1x3 block second
            ----------------------
*/
            DVdot13(nrowU, rowL0, Utemp0, Utemp1, Utemp2, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base2 += ichv2 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base2[-loc0] -= sums[2] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base2 -= ichv2 ;
         }
         offset += 3*nrowU ;
      }
      if ( irowL == lastInU - 1 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d and %d",
                 colindU[irowL], colindU[irowL+1]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + offset ;
         rowL1 = rowL0 + nrowU  ;
         colU0 = entU  + offset ;
         colU1 = colU0 + nrowU  ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
         ichv1 = colindU[irowL+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) ;
/*
         -----------------------------------
                 [ Ltemp0 ]   [ rowL0 ]
         compute [ Ltemp1 ] = [ rowL1 ] * D
         -----------------------------------
*/
         DVzero(2*nrowU, Ltemp0) ;
         SubMtx_scale2vec(mtxD, Ltemp0, Ltemp1, rowL0, rowL1) ;
/*
         -----------------------------------
                 [ Utemp0 ]       [ colU0 ]
         compute [ Utemp1 ] = D * [ colU0 ]
         -----------------------------------
*/
         DVzero(2*nrowU, Utemp0) ;
         SubMtx_scale2vec(mtxD, Utemp0, Utemp1, colU0, colU1) ;
/*
         ------------------------------------------------------------
         update the 2x2 diagonal block for these two rows and columns
         ------------------------------------------------------------
*/
         DVdot22(nrowU, Ltemp0, Ltemp1, colU0, colU1, sums) ;
/*
         -------------------------------------
         inject the four sums into the chevron
         -------------------------------------
*/
         base0[0]    -= sums[ 0] ; 
         loc = ichv1 - ichv0 ;
         base0[loc]  -= sums[ 1] ; 
         base0[-loc] -= sums[ 2] ; 
         base1[0]    -= sums[ 3] ; 
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL1 + nrowU ;
         colU0 = colU1 + nrowU ;
         for ( jcolU = irowL + 2 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + nrowU ;
            rowL2 = rowL1 + nrowU ;
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
            loc2 = colindU[jcolU+2] ;
/*
            ---------------------
            upper 2x3 block first
            ---------------------
*/
            DVdot23(nrowU, Ltemp0, Ltemp1, colU0, colU1, colU2, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base0[loc2] -= sums[2] ; 
            base1[loc0] -= sums[3] ; 
            base1[loc1] -= sums[4] ; 
            base1[loc2] -= sums[5] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
/*
            ----------------------
            lower 3x2 block second
            ----------------------
*/
            DVdot32(nrowU, rowL0, rowL1, rowL2, Utemp0, Utemp1, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base0[-loc1] -= sums[2] ; 
            base1[-loc1] -= sums[3] ; 
            base0[-loc2] -= sums[4] ; 
            base1[-loc2] -= sums[5] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + nrowU ;
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + nrowU ;
            colU1 = colU0 + nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
/*
            ---------------------
            upper 2x2 block first
            ---------------------
*/
            DVdot22(nrowU, Ltemp0, Ltemp1, colU0, colU1, sums) ;
/*
            -------------------------------------
            inject the four sums into the chevron
            -------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base1[loc0] -= sums[2] ; 
            base1[loc1] -= sums[3] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
/*
            ----------------------
            lower 2x2 block second
            ----------------------
*/
            DVdot22(nrowU, rowL0, rowL1, Utemp0, Utemp1, sums) ;
/*
            -------------------------------------
            inject the four sums into the chevron
            -------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base0[-loc1] -= sums[2] ; 
            base1[-loc1] -= sums[3] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = colindU[jcolU] ;
/*
            ---------------------
            upper 2x1 block first
            ---------------------
*/
            DVdot21(nrowU, Ltemp0, Ltemp1, colU0, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 -= ichv0 ;
            base1 -= ichv1 ;
            base0[loc0] -= sums[0] ; 
            base1[loc0] -= sums[1] ; 
            base0 += ichv0 ;
            base1 += ichv1 ;
/*
            ----------------------
            lower 1x2 block second
            ----------------------
*/
            DVdot12(nrowU, rowL0, Utemp0, Utemp1, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 += ichv0 ;
            base1 += ichv1 ;
            base0[-loc0] -= sums[0] ; 
            base1[-loc0] -= sums[1] ; 
            base0 -= ichv0 ;
            base1 -= ichv1 ;
         }
      } else if ( irowL == lastInU ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on row %d", colindU[irowL]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + offset ;
         colU0 = entU  + offset ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
/*
         -----------------------------------
         compute [ Ltemp0 ] = [ rowL0 ] * D
         -----------------------------------
*/
         DVzero(nrowU, Ltemp0) ;
         SubMtx_scale1vec(mtxD, Ltemp0, rowL0) ;
/*
         -----------------------------------
         compute [ Utemp0 ] = D * [ colU0 ]
         -----------------------------------
*/
         DVzero(nrowU, Utemp0) ;
         SubMtx_scale1vec(mtxD, Utemp0, colU0) ;
/*
         ------------------------------------------------------------
         update the 1x1 diagonal block for these two rows and columns
         ------------------------------------------------------------
*/
         DVdot11(nrowU, Ltemp0, colU0, sums) ;
/*
         -------------------------------
         inject the sum into the chevron
         -------------------------------
*/
         base0[0] -= sums[0] ; 
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL0 + nrowU ;
         colU0 = colU0 + nrowU ;
         for ( jcolU = irowL + 1 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + nrowU ;
            rowL2 = rowL1 + nrowU ;
            colU1 = colU0 + nrowU ;
            colU2 = colU1 + nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
            loc2 = colindU[jcolU+2] ;
/*
            ---------------------
            upper 1x3 block first
            ---------------------
*/
            DVdot13(nrowU, Ltemp0, colU0, colU1, colU2, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 -= ichv0 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base0[loc2] -= sums[2] ; 
            base0 += ichv0 ;
/*
            ----------------------
            lower 3x1 block second
            ----------------------
*/
            DVdot31(nrowU, rowL0, rowL1, rowL2, Utemp0, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 += ichv0 ;
            base0[-loc0] -= sums[0] ; 
            base0[-loc1] -= sums[1] ; 
            base0[-loc2] -= sums[2] ; 
            base0 -= ichv0 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + nrowU ;
            colU0 = colU2 + nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + nrowU ;
            colU1 = colU0 + nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = colindU[jcolU] ;
            loc1 = colindU[jcolU+1] ;
/*
            ---------------------
            upper 1x2 block first
            ---------------------
*/
            DVdot12(nrowU, Ltemp0, colU0, colU1, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 -= ichv0 ;
            base0[loc0] -= sums[0] ; 
            base0[loc1] -= sums[1] ; 
            base0 += ichv0 ;
/*
            ----------------------
            lower 2x1 block second
            ----------------------
*/
            DVdot21(nrowU, rowL0, rowL1, Utemp0, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 += ichv0 ;
            base0[-loc0] -= sums[0] ; 
            base0[-loc1] -= sums[1] ; 
            base0 -= ichv0 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = colindU[jcolU] ;
/*
            ---------------------
            upper 1x1 block first
            ---------------------
*/
            DVdot11(nrowU, Ltemp0, colU0, sums) ;
/*
            -------------------------------
            inject the sum into the chevron
            -------------------------------
*/
            base0 -= ichv0 ;
            base0[loc0] -= sums[0] ; 
            base0 += ichv0 ;
/*
            ----------------------
            lower 1x1 block second
            ----------------------
*/
            DVdot11(nrowU, rowL0, Utemp0, sums) ;
/*
            -------------------------------
            inject the sum into the chevron
            -------------------------------
*/
            base0 += ichv0 ;
            base0[-loc0] -= sums[0] ; 
            base0 -= ichv0 ;
         }
      }
   } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxU) 
            && SUBMTX_IS_SPARSE_ROWS(mtxL) ) {
      double   *base, *colU0, *colU1, *entL, *entU, *rowL0, *rowL1,
               *temp0, *temp1, *temp2 ;
      int      ichv, ii, irow0, irow1, jj, loc, ncolL, ncolU, nentL, 
               nentU, nrowL, nrowU, offsetL, offsetU, sizeL0, sizeL1, 
               sizeU0, sizeU1 ;
      int      *indL, *indL0, *indL1, *indU, *indU0, *indU1, 
               *sizesL, *sizesU ;

      SubMtx_sparseColumnsInfo(mtxU, 
                               &ncolU, &nentU, &sizesU, &indU, &entU) ;
      SubMtx_sparseRowsInfo(mtxL, 
                            &nrowL, &nentL, &sizesL, &indL, &entL) ;
      nrowU = mtxU->nrow ;
      ncolL = mtxL->ncol ;
      DV_setSize(tempDV, 3*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + nrowU ;
      temp2 = temp1 + nrowU ;
/*
      --------------------------------------------
      get the offsets into the indices and entries
      --------------------------------------------
*/
      for ( jcolU = offsetL = offsetU = 0 ; 
            jcolU < firstInU ; 
            jcolU++ ) {
         offsetL += sizesL[jcolU] ;
         offsetU += sizesU[jcolU] ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n %% offsetL %d, offsetU %d", 
              offsetL, offsetU) ;
      fflush(stdout) ;
#endif
/*
      ---------------------------------------------------
      loop over the supporting rows of L and columns of U
      ---------------------------------------------------
*/
      for ( irow0 = firstInU ; irow0 <= lastInU ; irow0++ ) {
         rowL0  = entL + offsetL ;
         indL0  = indL + offsetL ;
         colU0  = entU + offsetU ;
         indU0  = indU + offsetU ;
         sizeL0 = sizesL[irow0] ;
         sizeU0 = sizesU[irow0] ;
#if MYDEBUG > 0
         fprintf(stdout, 
       "\n\n %% irow0 %d, offsetL %d, offsetU %d, sizeL0 %d, sizeU0 %d",
              irow0, offsetL, offsetU, sizeL0, sizeU0) ;
         fflush(stdout) ;
#endif
         if ( sizeL0 > 0 || sizeU0 > 0 ) {
/*
            -----------------------------------------------------
            get base locations for the chevron's diagonal entries
            -----------------------------------------------------
*/
            ichv = colindU[irow0] ;
            base = Chv_diagLocation(chvT, ichv) ;
#if MYDEBUG > 0
         fprintf(stdout, "\n\n %% ichv = %d, base = %p", ichv, base) ;
         fflush(stdout) ;
#endif
            if ( sizeL0 > 0 ) {
/*
               ----------------------------------
               compute [ temp0 ] = D * [ rowL0 ]
               ----------------------------------
*/
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% loading temp1") ;
               fflush(stdout) ;
#endif
               DVzero(2*nrowU, temp0) ;
               for ( ii = 0 ; ii < sizeL0 ; ii++ ) {
                  jj = indL0[ii] ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% ii = %d, jj = %d", ii, jj) ;
                  fflush(stdout) ;
#endif
                  temp1[jj] = rowL0[ii] ;
               }
               SubMtx_scale1vec(mtxD, temp0, temp1) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% temp0 = L * D") ;
               DVfprintf(stdout, 2*nrowU, temp0) ;
               fflush(stdout) ;
#endif
            }
            if ( sizeU0 > 0 ) {
/*
               ----------------------------------
               compute [ temp1 ] = D * [ colU0 ]
               ----------------------------------
*/
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% loading temp2") ;
               fflush(stdout) ;
#endif
               DVzero(2*nrowU, temp1) ;
               for ( ii = 0 ; ii < sizeU0 ; ii++ ) {
                  jj = indU0[ii] ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% ii = %d, jj = %d", ii, jj) ;
                  fflush(stdout) ;
#endif
                  temp2[jj] = colU0[ii] ;
               }
               SubMtx_scale1vec(mtxD, temp1, temp2) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% temp1 = D * U") ;
               DVfprintf(stdout, nrowU, temp1) ;
               fflush(stdout) ;
#endif
            }
            if ( sizeL0 > 0 && sizeU0 > 0 ) {
/*
               -------------------------
               update the diagonal entry
               -------------------------
*/
               base[0] -= DVdoti(sizeU0, temp0, indU0, colU0) ;
            }
/*
            ----------------------------------------
            loop over the following rows and columns
            ----------------------------------------
*/
            if ( sizeU0 > 0 ) {
/*
               ------------------------
               update the lower entries
               ------------------------
*/
               base += ichv ;
               rowL1 = rowL0 + sizeL0 ;
               indL1 = indL0 + sizeL0 ;
               for ( irow1 = irow0 + 1 ; irow1 < ncolU ; irow1++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% irow1 %d, sizeL1 %d",
                          irow1, sizesL[irow1]) ;
                  fflush(stdout) ;
#endif
                  if ( (sizeL1 = sizesL[irow1]) > 0 ) {
                     loc = colindU[irow1] ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                       "\n\n %% base[%d] = %12.4e, base[%d] = %12.4e",
                       -loc, base[-loc], -loc + 1, base[-loc+1]) ;
                     fflush(stdout) ;
#endif
                     base[-loc] -= DVdoti(sizeL1, temp1, indL1, rowL1) ;
                     rowL1 += sizeL1 ;
                     indL1 += sizeL1 ;
                  }
               }
               base -= ichv ;
            }
            if ( sizeL0 > 0 ) {
/*
               ------------------------
               update the upper entries
               ------------------------
*/
               base -= ichv ;
               colU1 = colU0 + sizeU0 ;
               indU1 = indU0 + sizeU0 ;
               for ( irow1 = irow0 + 1 ; irow1 < ncolU ; irow1++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% irow1 %d, sizeU1 %d",
                          irow1, sizesU[irow1]) ;
                  fflush(stdout) ;
#endif
                  if ( (sizeU1 = sizesU[irow1]) > 0 ) {
                     loc = colindU[irow1] ;
                     base[loc] -= DVdoti(sizeU1, temp0, indU1, colU1) ;
                     colU1 += sizeU1 ;
                     indU1 += sizeU1 ;
                  }
               }
               base -= ichv ;
            }
         }
         offsetL += sizeL0 ;
         offsetU += sizeU0 ;
      }
   } else {
      fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p)"
              "\n bad input, mtxU must have dense or sparse columns"
              "\n and mtxL must have dense or sparse rows\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chvT) ) {
   if ( SUBMTX_IS_DENSE_COLUMNS(mtxU) && SUBMTX_IS_DENSE_ROWS(mtxL) ) {
      double   sums[18] ;
      double   *base0, *base1, *base2, *colU0, *colU1, *colU2, *entL,
               *entU, *Ltemp0, *Ltemp1, *Ltemp2, *rowL0, *rowL1, *rowL2,
               *Utemp0, *Utemp1, *Utemp2 ;
      int      ichv0, ichv1, ichv2, inc1, inc2, irowL, jcolU, 
               loc, loc0, loc1, loc2, ncolL, nrowL, nrowU, offset ;

      SubMtx_denseInfo(mtxL, &nrowL, &ncolL, &inc1, &inc2, &entL) ;
      SubMtx_denseInfo(mtxU, &nrowU, &ncolU, &inc1, &inc2, &entU) ;
      DV_setSize(tempDV, 12*nrowU) ;
      Ltemp0 = DV_entries(tempDV) ;
      Ltemp1 = Ltemp0 + 2*nrowU ;
      Ltemp2 = Ltemp1 + 2*nrowU ;
      Utemp0 = Ltemp2 + 2*nrowU ;
      Utemp1 = Utemp0 + 2*nrowU ;
      Utemp2 = Utemp1 + 2*nrowU ;
/*
      -----------------------------------------------------
      loop over the rows of L in groups of three
      updating the diagonal and upper blocks of the chevron
      -----------------------------------------------------
*/
      offset = firstInU*nrowU ;
      for ( irowL = firstInU ; irowL <= lastInU - 2 ; irowL += 3 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d, %d and %d",
                 colindU[irowL], colindU[irowL+1], colindU[irowL+2]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + 2*offset ;
         rowL1 = rowL0 + 2*nrowU  ;
         rowL2 = rowL1 + 2*nrowU  ;
         colU0 = entU  + 2*offset ;
         colU1 = colU0 + 2*nrowU  ;
         colU2 = colU1 + 2*nrowU  ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
         ichv1 = colindU[irowL+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) ;
         ichv2 = colindU[irowL+2] ;
         base2 = Chv_diagLocation(chvT, ichv2) ;
/*
         -----------------------------------
                 [ Ltemp0 ]   [ rowL0 ]
         compute [ Ltemp1 ] = [ rowL1 ] * D
                 [ Ltemp2 ]   [ rowL2 ]
         -----------------------------------
*/
         DVzero(6*nrowU, Ltemp0) ;
         SubMtx_scale3vec(mtxD, Ltemp0, Ltemp1, Ltemp2,
                        rowL0, rowL1, rowL2) ;
/*
         -----------------------------------
                 [ Utemp0 ]       [ colU0 ]
         compute [ Utemp1 ] = D * [ colU0 ]
                 [ Utemp2 ]       [ colU0 ]
         -----------------------------------
*/
         DVzero(6*nrowU, Utemp0) ;
         SubMtx_scale3vec(mtxD, Utemp0, Utemp1, Utemp2,
                        colU0, colU1, colU2) ;
/*
         --------------------------------------------------------------
         update the 3x3 diagonal block for these three rows and columns
         --------------------------------------------------------------
*/
         ZVdotU33(nrowU, Ltemp0, Ltemp1, Ltemp2, 
                  colU0, colU1, colU2, sums) ;
/*
         -------------------------------------
         inject the nine sums into the chevron
         -------------------------------------
*/
         base0[0]    -= sums[ 0] ; base0[1]      -= sums[ 1] ;
         loc = 2*(ichv1 - ichv0) ;
         base0[loc]  -= sums[ 2] ; base0[loc+1]  -= sums[ 3] ;
         base0[-loc] -= sums[ 6] ; base0[-loc+1] -= sums[ 7] ;
         loc = 2*(ichv2 - ichv0) ;
         base0[loc]  -= sums[ 4] ; base0[loc+1]  -= sums[ 5] ;
         base0[-loc] -= sums[12] ; base0[-loc+1] -= sums[13] ;
         base1[0]    -= sums[ 8] ; base1[1]      -= sums[ 9] ;
         loc = 2*(ichv2 - ichv1) ;
         base1[loc]  -= sums[10] ; base1[loc+1]  -= sums[11] ;
         base1[-loc] -= sums[14] ; base1[-loc+1] -= sums[15] ;
         base2[0]    -= sums[16] ; base2[1]      -= sums[17] ;
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL2 + 2*nrowU ;
         colU0 = colU2 + 2*nrowU ;
         for ( jcolU = irowL + 3 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + 2*nrowU ;
            rowL2 = rowL1 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
            loc2 = 2*colindU[jcolU+2] ;
/*
            ---------------------
            upper 3x3 block first
            ---------------------
*/
            ZVdotU33(nrowU, Ltemp0, Ltemp1, Ltemp2, 
                     colU0, colU1, colU2, sums) ;
/*
            -------------------------------------
            inject the nine sums into the chevron
            -------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base0[loc2] -= sums[ 4] ; base0[loc2+1] -= sums[ 5] ;
            base1[loc0] -= sums[ 6] ; base1[loc0+1] -= sums[ 7] ;
            base1[loc1] -= sums[ 8] ; base1[loc1+1] -= sums[ 9] ;
            base1[loc2] -= sums[10] ; base1[loc2+1] -= sums[11] ;
            base2[loc0] -= sums[12] ; base2[loc0+1] -= sums[13] ;
            base2[loc1] -= sums[14] ; base2[loc1+1] -= sums[15] ;
            base2[loc2] -= sums[16] ; base2[loc2+1] -= sums[17] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
/*
            ----------------------
            lower 3x3 block second
            ----------------------
*/
            ZVdotU33(nrowU, rowL0, rowL1, rowL2, 
                     Utemp0, Utemp1, Utemp2, sums) ;
/*
            -------------------------------------
            inject the nine sums into the chevron
            -------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base2[-loc0] -= sums[ 4] ; base2[-loc0+1] -= sums[ 5] ;
            base0[-loc1] -= sums[ 6] ; base0[-loc1+1] -= sums[ 7] ;
            base1[-loc1] -= sums[ 8] ; base1[-loc1+1] -= sums[ 9] ;
            base2[-loc1] -= sums[10] ; base2[-loc1+1] -= sums[11] ;
            base0[-loc2] -= sums[12] ; base0[-loc2+1] -= sums[13] ;
            base1[-loc2] -= sums[14] ; base1[-loc2+1] -= sums[15] ;
            base2[-loc2] -= sums[16] ; base2[-loc2+1] -= sums[17] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + 2*nrowU ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
/*
            ---------------------
            upper 3x2 block first
            ---------------------
*/
            ZVdotU32(nrowU, Ltemp0, Ltemp1, Ltemp2, 
                     colU0, colU1, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base1[loc0] -= sums[ 4] ; base1[loc0+1] -= sums[ 5] ;
            base1[loc1] -= sums[ 6] ; base1[loc1+1] -= sums[ 7] ;
            base2[loc0] -= sums[ 8] ; base2[loc0+1] -= sums[ 9] ;
            base2[loc1] -= sums[10] ; base2[loc1+1] -= sums[11] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
/*
            ----------------------
            lower 2x3 block second
            ----------------------
*/
            ZVdotU23(nrowU, rowL0, rowL1, 
                     Utemp0, Utemp1, Utemp2, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base2[-loc0] -= sums[ 4] ; base2[-loc0+1] -= sums[ 5] ;
            base0[-loc1] -= sums[ 6] ; base0[-loc1+1] -= sums[ 7] ;
            base1[-loc1] -= sums[ 8] ; base1[-loc1+1] -= sums[ 9] ;
            base2[-loc1] -= sums[10] ; base2[-loc1+1] -= sums[11] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
/*
            ---------------------
            upper 3x1 block first
            ---------------------
*/
            ZVdotU31(nrowU, Ltemp0, Ltemp1, Ltemp2, colU0, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base1[loc0] -= sums[ 2] ; base1[loc0+1] -= sums[ 3] ;
            base2[loc0] -= sums[ 4] ; base2[loc0+1] -= sums[ 5] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
/*
            ----------------------
            lower 1x3 block second
            ----------------------
*/
            ZVdotU13(nrowU, rowL0, Utemp0, Utemp1, Utemp2, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base2 += 2*ichv2 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base2[-loc0] -= sums[ 4] ; base2[-loc0+1] -= sums[ 5] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base2 -= 2*ichv2 ;
         }
         offset += 3*nrowU ;
      }
      if ( irowL == lastInU - 1 ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on rows %d and %d",
                 colindU[irowL], colindU[irowL+1]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + 2*offset ;
         rowL1 = rowL0 + 2*nrowU  ;
         colU0 = entU  + 2*offset ;
         colU1 = colU0 + 2*nrowU  ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
         ichv1 = colindU[irowL+1] ;
         base1 = Chv_diagLocation(chvT, ichv1) ;
/*
         -----------------------------------
                 [ Ltemp0 ]   [ rowL0 ]
         compute [ Ltemp1 ] = [ rowL1 ] * D
         -----------------------------------
*/
         DVzero(4*nrowU, Ltemp0) ;
         SubMtx_scale2vec(mtxD, Ltemp0, Ltemp1, rowL0, rowL1) ;
/*
         -----------------------------------
                 [ Utemp0 ]       [ colU0 ]
         compute [ Utemp1 ] = D * [ colU0 ]
         -----------------------------------
*/
         DVzero(4*nrowU, Utemp0) ;
         SubMtx_scale2vec(mtxD, Utemp0, Utemp1, colU0, colU1) ;
/*
         ------------------------------------------------------------
         update the 2x2 diagonal block for these two rows and columns
         ------------------------------------------------------------
*/
         ZVdotU22(nrowU, Ltemp0, Ltemp1, colU0, colU1, sums) ;
/*
         -------------------------------------
         inject the four sums into the chevron
         -------------------------------------
*/
         base0[0]    -= sums[ 0] ; base0[1]      -= sums[ 1] ;
         loc = 2*(ichv1 - ichv0) ;
         base0[loc]  -= sums[ 2] ; base0[loc+1]  -= sums[ 3] ;
         base0[-loc] -= sums[ 4] ; base0[-loc+1] -= sums[ 5] ;
         base1[0]    -= sums[ 6] ; base1[1]      -= sums[ 7] ;
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL1 + 2*nrowU ;
         colU0 = colU1 + 2*nrowU ;
         for ( jcolU = irowL + 2 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + 2*nrowU ;
            rowL2 = rowL1 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
            loc2 = 2*colindU[jcolU+2] ;
/*
            ---------------------
            upper 2x3 block first
            ---------------------
*/
            ZVdotU23(nrowU, Ltemp0, Ltemp1, colU0, colU1, colU2, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base0[loc2] -= sums[ 4] ; base0[loc2+1] -= sums[ 5] ;
            base1[loc0] -= sums[ 6] ; base1[loc0+1] -= sums[ 7] ;
            base1[loc1] -= sums[ 8] ; base1[loc1+1] -= sums[ 9] ;
            base1[loc2] -= sums[10] ; base1[loc2+1] -= sums[11] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
/*
            ----------------------
            lower 3x2 block second
            ----------------------
*/
            ZVdotU32(nrowU, rowL0, rowL1, rowL2, Utemp0, Utemp1, sums) ;
/*
            ------------------------------------
            inject the six sums into the chevron
            ------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base0[-loc1] -= sums[ 4] ; base0[-loc1+1] -= sums[ 5] ;
            base1[-loc1] -= sums[ 6] ; base1[-loc1+1] -= sums[ 7] ;
            base0[-loc2] -= sums[ 8] ; base0[-loc2+1] -= sums[ 9] ;
            base1[-loc2] -= sums[10] ; base1[-loc2+1] -= sums[11] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + 2*nrowU ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
/*
            ---------------------
            upper 2x2 block first
            ---------------------
*/
            ZVdotU22(nrowU, Ltemp0, Ltemp1, colU0, colU1, sums) ;
/*
            -------------------------------------
            inject the four sums into the chevron
            -------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base1[loc0] -= sums[ 4] ; base1[loc0+1] -= sums[ 5] ;
            base1[loc1] -= sums[ 6] ; base1[loc1+1] -= sums[ 7] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
/*
            ----------------------
            lower 2x2 block second
            ----------------------
*/
            ZVdotU22(nrowU, rowL0, rowL1, Utemp0, Utemp1, sums) ;
/*
            -------------------------------------
            inject the four sums into the chevron
            -------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base0[-loc1] -= sums[ 4] ; base0[-loc1+1] -= sums[ 5] ;
            base1[-loc1] -= sums[ 6] ; base1[-loc1+1] -= sums[ 7] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
/*
            ---------------------
            upper 2x1 block first
            ---------------------
*/
            ZVdotU21(nrowU, Ltemp0, Ltemp1, colU0, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base1[loc0] -= sums[ 2] ; base1[loc0+1] -= sums[ 3] ;
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
/*
            ----------------------
            lower 1x2 block second
            ----------------------
*/
            ZVdotU12(nrowU, rowL0, Utemp0, Utemp1, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 += 2*ichv0 ;
            base1 += 2*ichv1 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base1[-loc0] -= sums[ 2] ; base1[-loc0+1] -= sums[ 3] ;
            base0 -= 2*ichv0 ;
            base1 -= 2*ichv1 ;
         }
      } else if ( irowL == lastInU ) {
#if MYDEBUG > 0
         fprintf(stdout, "\n %% working on row %d", colindU[irowL]) ;
         fflush(stdout) ;
#endif
         rowL0 = entL  + 2*offset ;
         colU0 = entU  + 2*offset ;
/*
         -----------------------------------------------------
         get base locations for the chevron's diagonal entries
         -----------------------------------------------------
*/
         ichv0 = colindU[irowL] ;
         base0 = Chv_diagLocation(chvT, ichv0) ;
/*
         -----------------------------------
         compute [ Ltemp0 ] = [ rowL0 ] * D
         -----------------------------------
*/
         DVzero(2*nrowU, Ltemp0) ;
         SubMtx_scale1vec(mtxD, Ltemp0, rowL0) ;
/*
         -----------------------------------
         compute [ Utemp0 ] = D * [ colU0 ]
         -----------------------------------
*/
         DVzero(2*nrowU, Utemp0) ;
         SubMtx_scale1vec(mtxD, Utemp0, colU0) ;
/*
         ------------------------------------------------------------
         update the 1x1 diagonal block for these two rows and columns
         ------------------------------------------------------------
*/
         ZVdotU11(nrowU, Ltemp0, colU0, sums) ;
/*
         -------------------------------
         inject the sum into the chevron
         -------------------------------
*/
         base0[0] -= sums[ 0] ; base0[1] -= sums[ 1] ;
/*
         ------------------------------------------
         update the lower and upper blocks together
         ------------------------------------------
*/
         rowL0 = rowL0 + 2*nrowU ;
         colU0 = colU0 + 2*nrowU ;
         for ( jcolU = irowL + 1 ; jcolU < ncolU - 2 ; jcolU += 3 ) {
            rowL1 = rowL0 + 2*nrowU ;
            rowL2 = rowL1 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
            colU2 = colU1 + 2*nrowU ;
/*
            ------------------------------------------------------
            get the local indices for these three rows and columns
            ------------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
            loc2 = 2*colindU[jcolU+2] ;
/*
            ---------------------
            upper 1x3 block first
            ---------------------
*/
            ZVdotU13(nrowU, Ltemp0, colU0, colU1, colU2, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 -= 2*ichv0 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base0[loc2] -= sums[ 4] ; base0[loc2+1] -= sums[ 5] ;
            base0 += 2*ichv0 ;
/*
            ----------------------
            lower 3x1 block second
            ----------------------
*/
            ZVdotU31(nrowU, rowL0, rowL1, rowL2, Utemp0, sums) ;
/*
            --------------------------------------
            inject the three sums into the chevron
            --------------------------------------
*/
            base0 += 2*ichv0 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base0[-loc1] -= sums[ 2] ; base0[-loc1+1] -= sums[ 3] ;
            base0[-loc2] -= sums[ 4] ; base0[-loc2+1] -= sums[ 5] ;
            base0 -= 2*ichv0 ;
/*
            -------------------------------------
            increment the row and column pointers
            -------------------------------------
*/
            rowL0 = rowL2 + 2*nrowU ;
            colU0 = colU2 + 2*nrowU ;
         }
         if ( jcolU == ncolU - 2 ) {
            rowL1 = rowL0 + 2*nrowU ;
            colU1 = colU0 + 2*nrowU ;
/*
            ----------------------------------------------------
            get the local indices for these two rows and columns
            ----------------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
            loc1 = 2*colindU[jcolU+1] ;
/*
            ---------------------
            upper 1x2 block first
            ---------------------
*/
            ZVdotU12(nrowU, Ltemp0, colU0, colU1, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 -= 2*ichv0 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0[loc1] -= sums[ 2] ; base0[loc1+1] -= sums[ 3] ;
            base0 += 2*ichv0 ;
/*
            ----------------------
            lower 2x1 block second
            ----------------------
*/
            ZVdotU21(nrowU, rowL0, rowL1, Utemp0, sums) ;
/*
            ------------------------------------
            inject the two sums into the chevron
            ------------------------------------
*/
            base0 += 2*ichv0 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base0[-loc1] -= sums[ 2] ; base0[-loc1+1] -= sums[ 3] ;
            base0 -= 2*ichv0 ;
         } else if ( jcolU == ncolU - 1 ) {
/*
            ---------------------------------------------
            get the local indices for this row and column
            ---------------------------------------------
*/
            loc0 = 2*colindU[jcolU] ;
/*
            ---------------------
            upper 1x1 block first
            ---------------------
*/
            ZVdotU11(nrowU, Ltemp0, colU0, sums) ;
/*
            -------------------------------
            inject the sum into the chevron
            -------------------------------
*/
            base0 -= 2*ichv0 ;
            base0[loc0] -= sums[ 0] ; base0[loc0+1] -= sums[ 1] ;
            base0 += 2*ichv0 ;
/*
            ----------------------
            lower 1x1 block second
            ----------------------
*/
            ZVdotU11(nrowU, rowL0, Utemp0, sums) ;
/*
            -------------------------------
            inject the sum into the chevron
            -------------------------------
*/
            base0 += 2*ichv0 ;
            base0[-loc0] -= sums[ 0] ; base0[-loc0+1] -= sums[ 1] ;
            base0 -= 2*ichv0 ;
         }
      }
   } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxU) 
            && SUBMTX_IS_SPARSE_ROWS(mtxL) ) {
      double   imag, real ;
      double   *base, *colU0, *colU1, *entL, *entU, *rowL0, *rowL1,
               *temp0, *temp1, *temp2 ;
      int      ichv, ii, irow0, irow1, jj, loc, ncolL, ncolU, nentL, 
               nentU, nrowL, nrowU, offsetL, offsetU, sizeL0, sizeL1, 
               sizeU0, sizeU1 ;
      int      *indL, *indL0, *indL1, *indU, *indU0, *indU1, 
               *sizesL, *sizesU ;

      SubMtx_sparseColumnsInfo(mtxU, 
                               &ncolU, &nentU, &sizesU, &indU, &entU) ;
      SubMtx_sparseRowsInfo(mtxL, 
                            &nrowL, &nentL, &sizesL, &indL, &entL) ;
      nrowU = mtxU->nrow ;
      ncolL = mtxL->ncol ;
      DV_setSize(tempDV, 6*nrowU) ;
      temp0 = DV_entries(tempDV) ;
      temp1 = temp0 + 2*nrowU ;
      temp2 = temp1 + 2*nrowU ;
/*
      --------------------------------------------
      get the offsets into the indices and entries
      --------------------------------------------
*/
      for ( jcolU = offsetL = offsetU = 0 ; 
            jcolU < firstInU ; 
            jcolU++ ) {
         offsetL += sizesL[jcolU] ;
         offsetU += sizesU[jcolU] ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n %% offsetL %d, offsetU %d", 
              offsetL, offsetU) ;
      fflush(stdout) ;
#endif
/*
      ---------------------------------------------------
      loop over the supporting rows of L and columns of U
      ---------------------------------------------------
*/
      for ( irow0 = firstInU ; irow0 <= lastInU ; irow0++ ) {
         rowL0  = entL + 2*offsetL ;
         indL0  = indL +   offsetL ;
         colU0  = entU + 2*offsetU ;
         indU0  = indU +   offsetU ;
         sizeL0 = sizesL[irow0] ;
         sizeU0 = sizesU[irow0] ;
#if MYDEBUG > 0
         fprintf(stdout, 
       "\n\n %% irow0 %d, offsetL %d, offsetU %d, sizeL0 %d, sizeU0 %d",
              irow0, offsetL, offsetU, sizeL0, sizeU0) ;
         fflush(stdout) ;
#endif
         if ( sizeL0 > 0 || sizeU0 > 0 ) {
/*
            -----------------------------------------------------
            get base locations for the chevron's diagonal entries
            -----------------------------------------------------
*/
            ichv = colindU[irow0] ;
            base = Chv_diagLocation(chvT, ichv) ;
#if MYDEBUG > 0
         fprintf(stdout, "\n\n %% ichv = %d, base = %p", ichv, base) ;
         fflush(stdout) ;
#endif
            if ( sizeL0 > 0 ) {
/*
               ----------------------------------
               compute [ temp0 ] = D * [ rowL0 ]
               ----------------------------------
*/
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% loading temp1") ;
               fflush(stdout) ;
#endif
               DVzero(4*nrowU, temp0) ;
               for ( ii = 0 ; ii < sizeL0 ; ii++ ) {
                  jj = indL0[ii] ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% ii = %d, jj = %d", ii, jj) ;
                  fflush(stdout) ;
#endif
                  temp1[2*jj]   = rowL0[2*ii] ;
                  temp1[2*jj+1] = rowL0[2*ii+1] ;
#if MYDEBUG > 0
                  fprintf(stdout, ", (%12.4e,%12.4e)",
                          rowL0[2*ii], rowL0[2*ii+1]) ;
                  fflush(stdout) ;
#endif
               }
               SubMtx_scale1vec(mtxD, temp0, temp1) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% temp0 = L * D") ;
               DVfprintf(stdout, 2*nrowU, temp0) ;
               fflush(stdout) ;
#endif
            }
            if ( sizeU0 > 0 ) {
/*
               ----------------------------------
               compute [ temp1 ] = D * [ colU0 ]
               ----------------------------------
*/
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% loading temp2") ;
               fflush(stdout) ;
#endif
               DVzero(4*nrowU, temp1) ;
               for ( ii = 0 ; ii < sizeU0 ; ii++ ) {
                  jj = indU0[ii] ;
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% ii = %d, jj = %d", ii, jj) ;
                  fflush(stdout) ;
#endif
                  temp2[2*jj]   = colU0[2*ii] ;
                  temp2[2*jj+1] = colU0[2*ii+1] ;
#if MYDEBUG > 0
                  fprintf(stdout, ", (%12.4e,%12.4e)",
                          colU0[2*ii], colU0[2*ii+1]) ;
                  fflush(stdout) ;
#endif
               }
               SubMtx_scale1vec(mtxD, temp1, temp2) ;
#if MYDEBUG > 0
               fprintf(stdout, "\n\n %% temp1 = D * U") ;
               DVfprintf(stdout, 2*nrowU, temp1) ;
               fflush(stdout) ;
#endif
            }
            if ( sizeL0 > 0 && sizeU0 > 0 ) {
/*
               -------------------------
               update the diagonal entry
               -------------------------
*/
               ZVdotiU(sizeU0, temp0, indU0, colU0, &real, &imag) ;
               base[0] -= real ; base[1] -= imag ;
#if MYDEBUG > 0
               fprintf(stdout, 
                       "\n\n %% sum = %12.4e + %12.4e*i, "
                       "\n base[%d] = %12.4e, base[%d] = %12.4e ",
                       real, imag, 0, base[0], 1, base[1]) ;
               fflush(stdout) ;
#endif
            }
/*
            ----------------------------------------
            loop over the following rows and columns
            ----------------------------------------
*/
            if ( sizeU0 > 0 ) {
/*
               ------------------------
               update the lower entries
               ------------------------
*/
               base += 2*ichv ;
               rowL1 = rowL0 + 2*sizeL0 ;
               indL1 = indL0 +   sizeL0 ;
               for ( irow1 = irow0 + 1 ; irow1 < ncolU ; irow1++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% irow1 %d, sizeL1 %d",
                          irow1, sizesL[irow1]) ;
                  fflush(stdout) ;
#endif
                  if ( (sizeL1 = sizesL[irow1]) > 0 ) {
                     loc = 2*colindU[irow1] ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                       "\n\n %% base[%d] = %12.4e, base[%d] = %12.4e",
                       -loc, base[-loc], -loc + 1, base[-loc+1]) ;
                     fflush(stdout) ;
#endif
                     ZVdotiU(sizeL1, temp1, indL1, rowL1, &real, &imag);
                     base[-loc] -= real ; base[-loc+1] -= imag ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                    "\n\n %% sum = %12.4e + %12.4e*i, "
                       "\n base[%d] = %12.4e, base[%d] = %12.4e ",
                       real, imag, -loc, base[-loc], 
                       -loc + 1, base[-loc+1]) ;
                     fflush(stdout) ;
#endif
                     rowL1 += 2*sizeL1 ;
                     indL1 +=   sizeL1 ;
                  }
               }
               base -= 2*ichv ;
            }
            if ( sizeL0 > 0 ) {
/*
               ------------------------
               update the upper entries
               ------------------------
*/
               base -= 2*ichv ;
               colU1 = colU0 + 2*sizeU0 ;
               indU1 = indU0 +   sizeU0 ;
               for ( irow1 = irow0 + 1 ; irow1 < ncolU ; irow1++ ) {
#if MYDEBUG > 0
                  fprintf(stdout, "\n\n %% irow1 %d, sizeU1 %d",
                          irow1, sizesU[irow1]) ;
                  fflush(stdout) ;
#endif
                  if ( (sizeU1 = sizesU[irow1]) > 0 ) {
                     loc = 2*colindU[irow1] ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                       "\n\n %% base[%d] = %12.4e, base[%d] = %12.4e",
                       loc, base[loc], loc + 1, base[loc+1]) ;
                     fflush(stdout) ;
#endif
                     ZVdotiU(sizeU1, temp0, indU1, colU1, &real, &imag);
                     base[loc] -= real ; base[loc+1] -= imag ;
#if MYDEBUG > 0
                     fprintf(stdout, 
                       "\n\n %% sum = %12.4e + %12.4e*i, "
                       "\n base[%d] = %12.4e, base[%d] = %12.4e ",
                    real, imag, loc, base[loc], loc + 1, base[loc+1]) ;
                     fflush(stdout) ;
#endif
                     colU1 += 2*sizeU1 ;
                     indU1 +=   sizeU1 ;
                  }
               }
               base -= ichv ;
            }
         }
         offsetL += sizeL0 ;
         offsetU += sizeU0 ;
      }
   } else {
      fprintf(stderr, "\n fatal error in Chv_updateN(%p,%p,%p,%p)"
              "\n bad input, mtxU must have dense or sparse columns"
              "\n and mtxL must have dense or sparse rows\n",
              chvT, mtxD, mtxU, tempDV) ;
      exit(-1) ;
   }
}
/*
   ---------------------------------------------------
   overwrite the local indices with the global indices
   ---------------------------------------------------
*/
for ( jcolU = firstInU ; jcolU < ncolU ; jcolU++ ) {
   colindU[jcolU] = colindT[colindU[jcolU]] ;
}
return ; }

/*--------------------------------------------------------------------*/
