/*  split.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- for each U_{J,bnd{J}} matrix, remove from hash table,
              split into their U_{J,K} submatrices and insert 
              into the hash table.

   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_splitUpperMatrices (
   FrontMtx   *frontmtx,
   int        msglvl,
   FILE       *msgFile
) {
SubMtx          *mtxUJ, *mtxUJJ, *mtxUJK ;
SubMtxManager   *manager ;
double          *entUJ, *entUJK ;
int             count, first, ii, inc1, inc2, jcol, jj, J, K, nbytes,
                ncolJ, ncolUJ, ncolUJK, nentUJ, nentUJK, neqns, nfront, 
                nJ, nrowUJ, nrowUJK, offset, v ;
int             *colindJ, *colindUJ, *colindUJK, *colmap, *indicesUJ,
                *indicesUJK, *locmap, *rowindUJ, *rowindUJK, *sizesUJ, 
                *sizesUJK ;
I2Ohash         *upperhash ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_splitUpperMatrices(%p,%d,%p)"
           "\n bad input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
nfront    = FrontMtx_nfront(frontmtx) ;
neqns     = FrontMtx_neqns(frontmtx) ;
upperhash = frontmtx->upperhash ;
manager   = frontmtx->manager   ;
/*
   -----------------------------------
   construct the column and local maps
   -----------------------------------
*/
colmap = IVinit(neqns, -1) ;
locmap = IVinit(neqns, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
      if ( ncolJ > 0 && colindJ != NULL ) {
         for ( ii = 0 ; ii < nJ ; ii++ ) {
            v = colindJ[ii] ;
            colmap[v] = J ;
            locmap[v] = ii ;
         } 
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n colmap[]") ;
   IVfprintf(msgFile, neqns, colmap) ;
   fprintf(msgFile, "\n\n locmap[]") ;
   IVfprintf(msgFile, neqns, locmap) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   move the U_{J,J} matrices into the hash table
   ---------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (mtxUJJ = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
      I2Ohash_insert(frontmtx->upperhash, J, J, mtxUJJ) ;
   }
}
/*
   ------------------------------------------------------------
   now split the U_{J,bnd{J}} matrices into U_{J,K} matrices.
   note: columns of U_{J,bnd{J}} are assumed to be in ascending
   order with respect to the column ordering of the matrix.
   ------------------------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   mtxUJ = FrontMtx_upperMtx(frontmtx, J, nfront) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n ### J = %d, mtxUJ = %p", J, mtxUJ) ;
      fflush(msgFile) ;
   }
   if ( mtxUJ != NULL ) {
      if ( msglvl > 2 ) {
         SubMtx_writeForHumanEye(mtxUJ, msgFile) ;
         fflush(msgFile) ;
      }
      SubMtx_columnIndices(mtxUJ, &ncolUJ, &colindUJ) ;
      SubMtx_rowIndices(mtxUJ, &nrowUJ, &rowindUJ) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n  column indices for J") ;
         IVfprintf(msgFile, ncolUJ, colindUJ) ;
         fprintf(msgFile, "\n  row indices for UJ") ;
         IVfprintf(msgFile, nrowUJ, rowindUJ) ;
         fflush(msgFile) ;
      }
      if ( (K = colmap[colindUJ[0]]) == colmap[colindUJ[ncolUJ-1]] ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n  front %d supports only %d", J, K) ;
            fflush(msgFile) ;
         }
/*
         -------------------------------------------------
         U_{J,bnd{J}} is one submatrix, bnd{J} \subseteq K
         set row and column indices and change column id
         -------------------------------------------------
*/
         IVramp(nrowUJ, rowindUJ, 0, 1) ;
         for ( ii = 0 ; ii < ncolUJ ; ii++ ) {
            colindUJ[ii] = locmap[colindUJ[ii]] ;
         }
         SubMtx_setFields(mtxUJ, mtxUJ->type, mtxUJ->mode, J, K,
                          mtxUJ->nrow, mtxUJ->ncol, mtxUJ->nent) ;
/*
         mtxUJ->colid = K ;
*/
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n ##  inserting U(%d,%d) ", J, K) ;
            SubMtx_writeForHumanEye(mtxUJ, msgFile) ;
            fflush(msgFile) ;
         }
         I2Ohash_insert(upperhash, J, K, (void *) mtxUJ) ;
      } else {
/*
         -----------------------------------
         split U_{J,bnd{J}} into submatrices
         -----------------------------------
*/
         nJ = FrontMtx_frontSize(frontmtx, J) ;
         if ( SUBMTX_IS_DENSE_COLUMNS(mtxUJ) ) {
            SubMtx_denseInfo(mtxUJ, 
                           &nrowUJ, &ncolUJ, &inc1, &inc2, &entUJ) ;
         } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxUJ) ) {
            SubMtx_sparseColumnsInfo(mtxUJ, &ncolUJ, &nentUJ, 
                                   &sizesUJ, &indicesUJ, &entUJ) ;
            offset = 0 ;
            count  = sizesUJ[0] ;
         }
         first = 0 ;
         K = colmap[colindUJ[0]] ;
         for ( jcol = 1 ; jcol <= ncolUJ ; jcol++ ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n jcol = %d", jcol) ;
               if ( jcol < ncolUJ ) {
                  fprintf(msgFile, ", colmap[%d] = %d", 
                          colindUJ[jcol], colmap[colindUJ[jcol]]);
               }
               fflush(msgFile) ;
            }
            if ( jcol == ncolUJ || K != colmap[colindUJ[jcol]] ) {
               ncolUJK = jcol - first ;
               if ( SUBMTX_IS_DENSE_COLUMNS(mtxUJ) ) {
                  nentUJK = nJ*ncolUJK ;
               } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxUJ) ) {
                  if ( count == 0 ) {
                     goto no_entries ;
                  }
                  nentUJK = count ;
               }
               nbytes = SubMtx_nbytesNeeded(mtxUJ->type, mtxUJ->mode,
                                            nJ, ncolUJK, nentUJK) ;
               if ( msglvl > 2 ) {
                  fprintf(msgFile, 
                          "\n ncolUJK %d, nentUJK %d, nbytes %d",
                          ncolUJK, nentUJK, nbytes) ;
                  fflush(msgFile) ;
               }
               mtxUJK = SubMtxManager_newObjectOfSizeNbytes(manager, 
                                                          nbytes) ;
               SubMtx_init(mtxUJK, mtxUJ->type, mtxUJ->mode, J, K,
                         nJ, ncolUJK, nentUJK) ;
               if ( SUBMTX_IS_DENSE_COLUMNS(mtxUJ) ) {
                  SubMtx_denseInfo(mtxUJK, 
                         &nrowUJK, &ncolUJK, &inc1, &inc2, &entUJK) ;
                  if ( FRONTMTX_IS_REAL(frontmtx) ) {
                     DVcopy(nentUJK, entUJK, entUJ + first*nJ) ;
                  } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
                     DVcopy(2*nentUJK, entUJK, entUJ + 2*first*nJ) ;
                  }
               } else if ( SUBMTX_IS_SPARSE_COLUMNS(mtxUJ) ) {
                  SubMtx_sparseColumnsInfo(mtxUJK, &ncolUJK, &nentUJK, 
                                   &sizesUJK, &indicesUJK, &entUJK) ;
                  IVcopy(ncolUJK, sizesUJK, sizesUJ + first) ;
                  IVcopy(nentUJK, indicesUJK, indicesUJ + offset) ;
                  if ( FRONTMTX_IS_REAL(frontmtx) ) {
                     DVcopy(nentUJK, entUJK, entUJ + offset) ;
                  } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
                     DVcopy(2*nentUJK, entUJK, entUJ + 2*offset) ;
                  }
                  count  =  0 ;
                  offset += nentUJK ;
               }
/*
               -------------------------------------
               initialize the row and column indices
               -------------------------------------
*/
               if ( msglvl > 2 ) {
                  fprintf(msgFile, "\n setting row and column indices");
                  fflush(msgFile) ;
               }
               SubMtx_rowIndices(mtxUJK, &nrowUJK, &rowindUJK) ;
               IVramp(nJ, rowindUJK, 0, 1) ;
               SubMtx_columnIndices(mtxUJK, &ncolUJK, &colindUJK) ;
               for ( ii = 0, jj = first ; ii < ncolUJK ; ii++, jj++ ) {
                  colindUJK[ii] = locmap[colindUJ[jj]] ;
               }
/*
               ----------------------------------
               insert U_{J,K} into the hash table
               ----------------------------------
*/
               if ( msglvl > 2 ) {
                   fprintf(msgFile, 
                           "\n\n ##  inserting U(%d,%d) ", J, K) ;
                   SubMtx_writeForHumanEye(mtxUJK, msgFile) ;
                   fflush(msgFile) ;
               }
               I2Ohash_insert(upperhash, J, K, (void *) mtxUJK) ;
/*
               -----------------------------------
               we jump to here if there were no
               entries to be stored in the matrix.
               -----------------------------------
*/
   no_entries :
/*
               ----------------------------------------------------
               reset first and K to new first location and front id
               ----------------------------------------------------
*/
               first = jcol ;
               if ( jcol < ncolUJ ) {
                  K = colmap[colindUJ[jcol]] ;
               }
            } 
            if ( jcol < ncolUJ && SUBMTX_IS_SPARSE_COLUMNS(mtxUJ) ) {
               count += sizesUJ[jcol] ;
            }
         }
/*
         --------------------------------------------
         give U_{J,bnd{J}} back to the matrix manager
         --------------------------------------------
*/
         SubMtxManager_releaseObject(manager, mtxUJ) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(colmap) ;
IVfree(locmap) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- for each L_{bnd{J},J} matrix, remove from hash table,
              split into their L_{K,J} submatrices and insert 
              into the hash table.

   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_splitLowerMatrices (
   FrontMtx   *frontmtx,
   int         msglvl,
   FILE        *msgFile
) {
SubMtx          *mtxLJ, *mtxLJJ, *mtxLKJ ;
SubMtxManager   *manager ;
double        *entLJ, *entLKJ ;
int           count, first, ii, inc1, inc2, irow, jj, J, K, nbytes,
              ncolLJ, ncolLKJ, nentLJ, nentLKJ, neqns, nfront, nJ, 
              nrowJ, nrowLJ, nrowLKJ, offset, v ;
int           *colindLJ, *colindLKJ, *rowmap, *indicesLJ, *indicesLKJ, 
              *locmap, *rowindJ, *rowindLJ, *rowindLKJ, *sizesLJ, 
              *sizesLKJ ;
I2Ohash       *lowerhash ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_splitLowerMatrices(%p,%d,%p)"
           "\n bad input\n", frontmtx, msglvl, msgFile) ;
   exit(-1) ;
}
nfront    = FrontMtx_nfront(frontmtx) ;
neqns     = FrontMtx_neqns(frontmtx) ;
lowerhash = frontmtx->lowerhash ;
manager   = frontmtx->manager   ;
/*
   --------------------------------
   construct the row and local maps
   --------------------------------
*/
rowmap = IVinit(neqns, -1) ;
locmap = IVinit(neqns, -1) ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (nJ = FrontMtx_frontSize(frontmtx, J)) > 0 ) {
      FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
      if ( nrowJ > 0 && rowindJ != NULL ) {
         for ( ii = 0 ; ii < nJ ; ii++ ) {
            v = rowindJ[ii] ;
            rowmap[v] = J ;
            locmap[v] = ii ;
         } 
      }
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rowmap[]") ;
   IVfprintf(msgFile, neqns, rowmap) ;
   fprintf(msgFile, "\n\n locmap[]") ;
   IVfprintf(msgFile, neqns, locmap) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   move the L_{J,J} matrices into the hash table
   ---------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   if ( (mtxLJJ = FrontMtx_lowerMtx(frontmtx, J, J)) != NULL ) {
      I2Ohash_insert(frontmtx->lowerhash, J, J, mtxLJJ) ;
   }
}
/*
   ------------------------------------------------------------
   now split the L_{bnd{J},J} matrices into L_{K,J} matrices.
   note: columns of L_{bnd{J},J} are assumed to be in ascending
   order with respect to the column ordering of the matrix.
   ------------------------------------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   mtxLJ = FrontMtx_lowerMtx(frontmtx, nfront, J) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n ### J = %d, mtxLJ = %p", J, mtxLJ) ;
      fflush(msgFile) ;
   }
   if ( mtxLJ != NULL ) {
      if ( msglvl > 2 ) {
         SubMtx_writeForHumanEye(mtxLJ, msgFile) ;
         fflush(msgFile) ;
      }
      SubMtx_columnIndices(mtxLJ, &ncolLJ, &colindLJ) ;
      SubMtx_rowIndices(mtxLJ, &nrowLJ, &rowindLJ) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n  column indices for J") ;
         IVfprintf(msgFile, ncolLJ, colindLJ) ;
         fprintf(msgFile, "\n  row indices for LJ") ;
         IVfprintf(msgFile, nrowLJ, rowindLJ) ;
         fflush(msgFile) ;
      }
      if ( (K = rowmap[rowindLJ[0]]) == rowmap[rowindLJ[nrowLJ-1]] ) {
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n  front %d supports only %d", J, K) ;
            fflush(msgFile) ;
         }
/*
         -------------------------------------------------
         L_{bnd{J},J} is one submatrix, bnd{J} \subseteq K
         set row and column indices and change column id
         -------------------------------------------------
*/
         IVramp(ncolLJ, colindLJ, 0, 1) ;
         for ( ii = 0 ; ii < nrowLJ ; ii++ ) {
            rowindLJ[ii] = locmap[rowindLJ[ii]] ;
         }
/*
         mtxLJ->rowid = K ;
*/
         SubMtx_setFields(mtxLJ, mtxLJ->type, mtxLJ->mode, K, J,
                          mtxLJ->nrow, mtxLJ->ncol, mtxLJ->nent) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n ##  inserting L(%d,%d) ", K, J) ;
            SubMtx_writeForHumanEye(mtxLJ, msgFile) ;
            fflush(msgFile) ;
         }
         I2Ohash_insert(lowerhash, K, J, (void *) mtxLJ) ;
      } else {
/*
         -----------------------------------
         split L_{bnd{J},J} into submatrices
         -----------------------------------
*/
         nJ = FrontMtx_frontSize(frontmtx, J) ;
         if ( SUBMTX_IS_DENSE_ROWS(mtxLJ) ) {
            SubMtx_denseInfo(mtxLJ, 
                           &nrowLJ, &ncolLJ, &inc1, &inc2, &entLJ) ;
         } else if ( SUBMTX_IS_SPARSE_ROWS(mtxLJ) ) {
            SubMtx_sparseRowsInfo(mtxLJ, &nrowLJ, &nentLJ, 
                                &sizesLJ, &indicesLJ, &entLJ) ;
            offset = 0 ;
            count  = sizesLJ[0] ;
         }
         first = 0 ;
         K = rowmap[rowindLJ[0]] ;
         for ( irow = 1 ; irow <= nrowLJ ; irow++ ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n irow = %d", irow) ;
               if ( irow < nrowLJ ) {
                  fprintf(msgFile, ", rowmap[%d] = %d", 
                          rowindLJ[irow], rowmap[rowindLJ[irow]]);
               }
               fflush(msgFile) ;
            }
            if ( irow == nrowLJ || K != rowmap[rowindLJ[irow]] ) {
               nrowLKJ = irow - first ;
               if ( SUBMTX_IS_DENSE_ROWS(mtxLJ) ) {
                  nentLKJ = nJ*nrowLKJ ;
               } else if ( SUBMTX_IS_SPARSE_ROWS(mtxLJ) ) {
                  if ( count == 0 ) {
                     goto no_entries ;
                  }
                  nentLKJ = count ;
               }
               nbytes = SubMtx_nbytesNeeded(mtxLJ->type, mtxLJ->mode,
                                            nrowLKJ, nJ, nentLKJ) ;
               mtxLKJ = SubMtxManager_newObjectOfSizeNbytes(manager, 
                                                          nbytes) ;
               SubMtx_init(mtxLKJ, mtxLJ->type, mtxLJ->mode, K, J,
                         nrowLKJ, nJ, nentLKJ) ;
               if ( SUBMTX_IS_DENSE_ROWS(mtxLJ) ) {
                  SubMtx_denseInfo(mtxLKJ, 
                         &nrowLKJ, &ncolLKJ, &inc1, &inc2, &entLKJ) ;
                  if ( FRONTMTX_IS_REAL(frontmtx) ) {
                     DVcopy(nentLKJ, entLKJ, entLJ + first*nJ) ;
                  } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
                     DVcopy(2*nentLKJ, entLKJ, entLJ + 2*first*nJ) ;
                  }
               } else if ( SUBMTX_IS_SPARSE_ROWS(mtxLJ) ) {
                  SubMtx_sparseRowsInfo(mtxLKJ, &nrowLKJ, &nentLKJ, 
                                      &sizesLKJ, &indicesLKJ, &entLKJ) ;
                  IVcopy(nrowLKJ, sizesLKJ, sizesLJ + first) ;
                  IVcopy(nentLKJ, indicesLKJ, indicesLJ + offset) ;
                  if ( FRONTMTX_IS_REAL(frontmtx) ) {
                     DVcopy(nentLKJ, entLKJ, entLJ + offset) ;
                  } else if ( FRONTMTX_IS_COMPLEX(frontmtx) ) {
                     DVcopy(2*nentLKJ, entLKJ, entLJ + 2*offset) ;
                  }
                  count  =  0 ;
                  offset += nentLKJ ;
               }
/*
               -------------------------------------
               initialize the row and column indices
               -------------------------------------
*/
               SubMtx_rowIndices(mtxLKJ, &nrowLKJ, &rowindLKJ) ;
               for ( ii = 0, jj = first ; ii < nrowLKJ ; ii++, jj++ ) {
                  rowindLKJ[ii] = locmap[rowindLJ[jj]] ;
               }
               SubMtx_columnIndices(mtxLKJ, &ncolLKJ, &colindLKJ) ;
               IVramp(ncolLKJ, colindLKJ, 0, 1) ;
/*
               ----------------------------------
               insert L_{K,J} into the hash table
               ----------------------------------
*/
               if ( msglvl > 2 ) {
                   fprintf(msgFile, 
                           "\n\n ##  inserting L(%d,%d) ", K, J) ;
                   SubMtx_writeForHumanEye(mtxLKJ, msgFile) ;
                   fflush(msgFile) ;
               }
               I2Ohash_insert(lowerhash, K, J, (void *) mtxLKJ) ;
/*
               -----------------------------------
               we jump to here if there were no
               entries to be stored in the matrix.
               -----------------------------------
*/
   no_entries :
/*
               ----------------------------------------------------
               reset first and K to new first location and front id
               ----------------------------------------------------
*/
               first = irow ;
               if ( irow < nrowLJ ) {
                  K = rowmap[rowindLJ[irow]] ;
               }
            } 
            if ( irow < nrowLJ && SUBMTX_IS_SPARSE_ROWS(mtxLJ) ) {
               count += sizesLJ[irow] ;
            }
         }
/*
         --------------------------------------------
         give L_{bnd{J},J} back to the matrix manager
         --------------------------------------------
*/
         SubMtxManager_releaseObject(manager, mtxLJ) ;
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(rowmap) ;
IVfree(locmap) ;

return ; }

/*--------------------------------------------------------------------*/
