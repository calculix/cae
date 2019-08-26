/*  util.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- produce a map from each column 
              to the front that contains it

   created -- 98may04, cca
   -----------------------------------------
*/
IV *
FrontMtx_colmapIV ( 
   FrontMtx   *frontmtx
) {
int   ii, J, ncolJ, neqns, nfront, nJ ;
int   *colindJ, *colmap ;
IV    *colmapIV ;
/*
   -----------------------------------------
   get the map from columns to owning fronts
   -----------------------------------------
*/
neqns  = FrontMtx_neqns(frontmtx) ;
nfront = FrontMtx_nfront(frontmtx) ;
colmapIV = IV_new() ;
IV_init(colmapIV, neqns, NULL) ;
colmap = IV_entries(colmapIV) ;
IVfill(neqns, colmap, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
      if ( ncolJ > 0 && colindJ != NULL ) {
         for ( ii = 0 ; ii < nJ ; ii++ ) {
            colmap[colindJ[ii]] = J ;
         }
      }
   }
}
return(colmapIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- produce a map from each row to the front that contains it

   created -- 98may04, cca
   --------------------------------------------------------------------
*/
IV *
FrontMtx_rowmapIV ( 
   FrontMtx   *frontmtx
) {
int   ii, J, nrowJ, neqns, nfront, nJ ;
int   *rowindJ, *rowmap ;
IV    *rowmapIV ;
/*
   --------------------------------------
   get the map from rows to owning fronts
   --------------------------------------
*/
neqns  = FrontMtx_neqns(frontmtx) ;
nfront = FrontMtx_nfront(frontmtx) ;
rowmapIV = IV_new() ;
IV_init(rowmapIV, neqns, NULL) ;
rowmap = IV_entries(rowmapIV) ;
IVfill(neqns, rowmap, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
      if ( nrowJ > 0 && rowindJ != NULL ) {
         for ( ii = 0 ; ii < nJ ; ii++ ) {
            rowmap[rowindJ[ii]] = J ;
         }
      }
   }
}
return(rowmapIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   compute the inertia of a symmetric matrix

   fill *pnnegative with the number of negative eigenvalues of A
   fill *pnzero     with the number of zero eigenvalues of A
   fill *pnpositive with the number of positive eigenvalues of A

   created -- 98may04, cca
   -------------------------------------------------------------
*/
void
FrontMtx_inertia (
   FrontMtx   *frontmtx,
   int        *pnnegative,
   int        *pnzero,
   int        *pnpositive
) {
SubMtx     *mtx ;
double   arm, areal, bimag, breal, creal, mid, val ;
double   *entries ;
int      ii, ipivot, irow, J, nent, nfront, nJ, 
         nnegative, npositive, nzero ;
int      *pivotsizes ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL 
   || pnnegative == NULL || pnzero == NULL || pnpositive == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_inertia(%p,%p,%p,%p)"
           "\n bad input\n", 
           frontmtx, pnnegative, pnzero, pnpositive) ;
   fflush(stdout) ;
}
if ( FRONTMTX_IS_REAL(frontmtx) && ! FRONTMTX_IS_SYMMETRIC(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_inertia(%p,%p,%p,%p)"
           "\n matrix is real and not symmetric \n",
           frontmtx, pnnegative, pnzero, pnpositive) ;
   fflush(stdout) ;
} else if ( FRONTMTX_IS_COMPLEX(frontmtx) 
        &&  ! FRONTMTX_IS_HERMITIAN(frontmtx) ) {
   fprintf(stderr, "\n fatal error in FrontMtx_inertia(%p,%p,%p,%p)"
           "\n matrix is complex and not hermitian \n",
           frontmtx, pnnegative, pnzero, pnpositive) ;
   fflush(stdout) ;
}
nfront = frontmtx->nfront ;
nnegative = nzero = npositive = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   mtx = FrontMtx_diagMtx(frontmtx, J) ;
   if ( mtx != NULL ) {
      if ( ! FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
         -----------
         no pivoting
         -----------
*/
         SubMtx_diagonalInfo(mtx, &nJ, &entries) ;
         if ( FRONTMTX_IS_REAL(frontmtx) ) {
            for ( ii = 0 ; ii < nJ ; ii++ ) {
               if ( entries[ii] < 0.0 ) {
                  nnegative++ ;
               } else if ( entries[ii] > 0.0 ) {
                  npositive++ ;
               } else {
                  nzero++ ;
               }
            }
         } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
            for ( ii = 0 ; ii < nJ ; ii++ ) {
               if ( entries[2*ii] < 0.0 ) {
                  nnegative++ ;
               } else if ( entries[2*ii] > 0.0 ) {
                  npositive++ ;
               } else {
                  nzero++ ;
               }
            }
         }
      } else {
/*
         --------
         pivoting
         --------
*/
         SubMtx_blockDiagonalInfo(mtx, &nJ, &nent, 
                                  &pivotsizes, &entries) ;
         if ( FRONTMTX_IS_REAL(frontmtx) ) {
            for ( irow = ipivot = ii = 0 ; irow < nJ ; ipivot++ ) {
               if ( pivotsizes[ipivot] == 1 ) {
                  val = entries[ii] ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  irow++ ; ii++ ;
               } else {
                  areal = entries[ii] ;
                  breal = entries[ii+1] ;
                  creal = entries[ii+2] ;
                  mid = 0.5*(areal + creal) ;
                  arm = sqrt(0.25*(areal - creal)*(areal - creal) 
                             + breal*breal) ;
                  val = mid + arm ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  val = mid - arm ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  irow += 2 ; ii += 3 ;
               }
            }
         } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
            for ( irow = ipivot = ii = 0 ; irow < nJ ; ipivot++ ) {
               if ( pivotsizes[ipivot] == 1 ) {
                  val = entries[2*ii] ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  irow++ ; ii++ ;
               } else {
                  areal = entries[2*ii]   ;
                  breal = entries[2*ii+2] ;
                  bimag = entries[2*ii+3] ;
                  creal = entries[2*ii+4] ;
                  mid = 0.5*(areal + creal) ;
                  arm = sqrt(0.25*(areal - creal)*(areal - creal) 
                             + breal*breal + bimag*bimag) ;
                  val = mid + arm ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  val = mid - arm ;
                  if ( val < 0.0 ) {
                     nnegative++ ;
                  } else if ( val > 0.0 ) {
                     npositive++ ;
                  } else {
                     nzero++ ;
                  }
                  irow += 2 ; ii += 3 ;
               }
            }
         }
      }
   }
}
*pnnegative = nnegative ;
*pnzero     = nzero     ;
*pnpositive = npositive ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- create and return an IV object that contains 
              all the row ids owned by process myid.

   created -- 98jun13, cca
   -------------------------------------------------------
*/
IV *
FrontMtx_ownedRowsIV (
   FrontMtx   *frontmtx,
   int        myid,
   IV         *ownersIV,
   int        msglvl,
   FILE       *msgFile
) {
int   J, neqns, nfront, nJ, nowned, nrowJ, offset ;
int   *ownedRows, *owners, *rowindJ ;
IV    *ownedRowsIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_ownedRowsIV(%p,%d,%p)"
           "\n bad input\n", frontmtx, myid, ownersIV) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
neqns  = frontmtx->neqns  ;
ownedRowsIV = IV_new() ;
if ( ownersIV == NULL ) {
   IV_init(ownedRowsIV, neqns, NULL) ;
   IV_ramp(ownedRowsIV, 0, 1) ;
} else {
   owners = IV_entries(ownersIV) ;
   for ( J = 0, nowned = 0 ; J < nfront ; J++ ) {
      if ( owners[J] == myid ) {
         nJ = FrontMtx_frontSize(frontmtx, J) ;
         nowned += nJ ;
      }
   }
   if ( nowned > 0 ) {
      IV_init(ownedRowsIV, nowned, NULL) ;
      ownedRows = IV_entries(ownedRowsIV) ;
      for ( J = 0, offset = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid ) {
            nJ = FrontMtx_frontSize(frontmtx, J) ;
            if ( nJ > 0 ) {
               FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
               IVcopy(nJ, ownedRows + offset, rowindJ) ;
               offset += nJ ;
            }
         }
      }
   }
}
return(ownedRowsIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- create and return an IV object that contains 
              all the column ids owned by process myid.

   created -- 98jun13, cca
   -------------------------------------------------------
*/
IV *
FrontMtx_ownedColumnsIV (
   FrontMtx   *frontmtx,
   int        myid,
   IV         *ownersIV,
   int        msglvl,
   FILE       *msgFile
) {
int   J, neqns, nfront, nJ, nowned, ncolJ, offset ;
int   *ownedColumns, *owners, *colindJ ;
IV    *ownedColumnsIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_ownedColumnsIV(%p,%d,%p)"
           "\n bad input\n", frontmtx, myid, ownersIV) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
neqns  = frontmtx->neqns  ;
ownedColumnsIV = IV_new() ;
if ( ownersIV == NULL ) {
   IV_init(ownedColumnsIV, neqns, NULL) ;
   IV_ramp(ownedColumnsIV, 0, 1) ;
} else {
   owners = IV_entries(ownersIV) ;
   for ( J = 0, nowned = 0 ; J < nfront ; J++ ) {
      if ( owners[J] == myid ) {
         nJ = FrontMtx_frontSize(frontmtx, J) ;
         nowned += nJ ;
      }
   }
   if ( nowned > 0 ) {
      IV_init(ownedColumnsIV, nowned, NULL) ;
      ownedColumns = IV_entries(ownedColumnsIV) ;
      for ( J = 0, offset = 0 ; J < nfront ; J++ ) {
         if ( owners[J] == myid ) {
            nJ = FrontMtx_frontSize(frontmtx, J) ;
            if ( nJ > 0 ) {
               FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
               IVcopy(nJ, ownedColumns + offset, colindJ) ;
               offset += nJ ;
            }
         }
      }
   }
}
return(ownedColumnsIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to create and return an IVL object that holds the
      submatrix nonzero pattern for the upper triangular factor.

   NOTE: this method supercedes calling IVL_mapEntries() on
         the column adjacency structure. that gave problems when
         pivoting forced some fronts to have no eliminated columns.
         in some cases, solve aggregates were expected when none
         were forthcoming.

   created -- 98aug20, cca
   ----------------------------------------------------------------
*/
IVL *
FrontMtx_makeUpperBlockIVL (
   FrontMtx   *frontmtx,
   IV         *colmapIV
) {
int   count, ii, J, K, ncol, nfront, nJ ;
int   *colmap, *colind, *list, *mark ;
IVL   *upperblockIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || colmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_makeUpperBlockIVL()"
           "\n bad input\n") ;
   exit(-1) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
colmap = IV_entries(colmapIV) ;
/*
   -----------------------------
   set up the working storage
   and initialize the IVL object
   -----------------------------
*/
mark = IVinit(nfront, -1) ;
list = IVinit(nfront, -1) ;
upperblockIVL = IVL_new() ;
IVL_init1(upperblockIVL, IVL_CHUNKED, nfront) ;
/*
   -------------------
   fill the IVL object
   -------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_columnIndices(frontmtx, J, &ncol, &colind) ;
      if ( ncol > 0 ) {
         mark[J] = J ;
         count = 0 ;
         list[count++] = J ;
         for ( ii = nJ ; ii < ncol ; ii++ ) {
            K = colmap[colind[ii]] ;
            if ( mark[K] != J ) {
               mark[K] = J ;
               list[count++] = K ;
            }
         }
         IVL_setList(upperblockIVL, J, count, list) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(list) ;

return(upperblockIVL) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to create and return an IVL object that holds the
      submatrix nonzero pattern for the lower triangular factor.

   NOTE: this method supercedes calling IVL_mapEntries() on
         the row adjacency structure. that gave problems when
         pivoting forced some fronts to have no eliminated columns.
         in some cases, solve aggregates were expected when none
         were forthcoming.

   created -- 98aug20, cca
   ----------------------------------------------------------------
*/
IVL *
FrontMtx_makeLowerBlockIVL (
   FrontMtx   *frontmtx,
   IV         *rowmapIV
) {
int   count, ii, J, K, nrow, nfront, nJ ;
int   *rowmap, *rowind, *list, *mark ;
IVL   *lowerblockIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || rowmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_makeLowerBlockIVL()"
           "\n bad input\n") ;
   exit(-1) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
rowmap = IV_entries(rowmapIV) ;
/*
   -----------------------------
   set up the working storage
   and initialize the IVL object
   -----------------------------
*/
mark = IVinit(nfront, -1) ;
list = IVinit(nfront, -1) ;
lowerblockIVL = IVL_new() ;
IVL_init1(lowerblockIVL, IVL_CHUNKED, nfront) ;
/*
   -------------------
   fill the IVL object
   -------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_rowIndices(frontmtx, J, &nrow, &rowind) ;
      if ( nrow > 0 ) {
         mark[J] = J ;
         count = 0 ;
         list[count++] = J ;
         for ( ii = nJ ; ii < nrow ; ii++ ) {
            K = rowmap[rowind[ii]] ;
            if ( mark[K] != J ) {
               mark[K] = J ;
               list[count++] = K ;
            }
         }
         IVL_setList(lowerblockIVL, J, count, list) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(list) ;

return(lowerblockIVL) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to compute and return the number of floating point
      operations to perform a solve with a single right hand side.

   created -- 98sep26, cca
   ---------------------------------------------------------------
*/
int
FrontMtx_nSolveOps (
   FrontMtx   *frontmtx
) {
int   nsolveops ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_nSolveOps()"
           "\n frontmtx is NULL\n") ;
   exit(-1) ;
}
switch ( frontmtx->type ) {
case SPOOLES_REAL :
   switch ( frontmtx->symmetryflag ) {
      case SPOOLES_SYMMETRIC :
         nsolveops = 4*frontmtx->nentU + frontmtx->nentD ;
         break ;
      case SPOOLES_NONSYMMETRIC :
         nsolveops = 2*frontmtx->nentL + frontmtx->nentD 
                      + 2*frontmtx->nentU ;
         break ;
      default :
         fprintf(stderr, "\n fatal error in FrontMtx_nSolveOps()"
                 "\n real type, invalid symmetryflag %d\n", 
                 frontmtx->symmetryflag) ;
         exit(-1) ;
         break ;
   }
   break ;
case SPOOLES_COMPLEX :
   switch ( frontmtx->symmetryflag ) {
      case SPOOLES_SYMMETRIC :
      case SPOOLES_HERMITIAN :
         nsolveops = 16*frontmtx->nentU + 8*frontmtx->nentD ;
         break ;
      case SPOOLES_NONSYMMETRIC :
         nsolveops = 8*frontmtx->nentL + 8*frontmtx->nentD 
                      + 8*frontmtx->nentU ;
         break ;
      default :
         fprintf(stderr, "\n fatal error in FrontMtx_nSolveOps()"
                 "\n complex type, invalid symmetryflag %d\n", 
                 frontmtx->symmetryflag) ;
         exit(-1) ;
         break ;
   }
   break ;
default :
   fprintf(stderr, "\n fatal error in FrontMtx_nSolveOps()"
           "\n invalid type %d\n", frontmtx->type) ;
   exit(-1) ;
   break ;
}
return(nsolveops) ; }

/*--------------------------------------------------------------------*/
