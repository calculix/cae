/*  storeFront.c  */

#include "../FrontMtx.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to store the factor indicies and entries from front J

   pivotsizesIV -- used for symmetric or hermitian and pivoting
   droptol      -- used for drop tolerance factorization,
                   an entry is stored if its magnitude > droptol

   created -- 98may04, cca
   ----------------------------------------------------------------
*/
void
FrontMtx_storeFront (
   FrontMtx   *frontmtx,
   Chv        *frontJ,
   IV         *pivotsizesIV,
   double     droptol,
   int        msglvl,
   FILE       *msgFile
) {
SubMtx          *mtx ;
double          *entries ;
int             inc1, inc2, ipivot, irow, jcol, J, m, nbytes, 
                ncol, nD, nent, nfront, nL, npivot, nrow, nU ;
int             *colind, *colindJ, *firstlocs, *indices, *pivots, 
                *pivotsizes, *rowind, *rowindJ, *sizes ;
SubMtxManager   *manager ;
/*
   -------------------------------------
   extract information from front matrix
   -------------------------------------
*/
nfront  = frontmtx->nfront  ;
manager = frontmtx->manager ;
if ( (   FRONTMTX_IS_SYMMETRIC(frontmtx) 
      || FRONTMTX_IS_HERMITIAN(frontmtx) )
   && FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IV_sizeAndEntries(pivotsizesIV, &npivot, &pivotsizes) ;
} else {
   npivot = 0, pivotsizes = NULL ;
}
/*
   --------------------------------
   extract information from chevron
   --------------------------------
*/
J = frontJ->id ;
Chv_dimensions(frontJ, &nD, &nL, &nU) ;
Chv_columnIndices(frontJ, &ncol, &colindJ) ;
Chv_rowIndices(frontJ, &nrow, &rowindJ) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n inside FrontMtx_storeFront, J = %d", J) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   if ( frontmtx->lock != NULL ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n locking lock") ;
         fflush(msgFile) ;
      }
      Lock_lock(frontmtx->lock) ;
      frontmtx->nlocks++ ;
   }
/*
   --------------------
   store the front size
   --------------------
*/
   FrontMtx_setFrontSize(frontmtx, J, nD) ;
/*
   ------------------------
   store the column indices
   ------------------------
*/
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n setting column indices, J = %d", J) ;
      fflush(msgFile) ;
   }
   IVL_setList(frontmtx->coladjIVL, J, ncol, colindJ) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      ---------------------
      store the row indices
      ---------------------
*/
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n setting row indices, J = %d", J) ;
         fflush(msgFile) ;
      }
      IVL_setList(frontmtx->rowadjIVL, J, nrow, rowindJ) ;
   }
   if ( frontmtx->lock != NULL ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n unlocking lock") ;
         fflush(msgFile) ;
      }
      Lock_unlock(frontmtx->lock) ;
   }
}
if ( nD == 0 ) {
   return ;
}
/*
   ------------------------
   store the D_{J,J} matrix
   ------------------------
*/
if ( pivotsizes != NULL ) {
/* 
   ---------------------------------------------------
   symmetric and pivoting, store block diagonal matrix
   ---------------------------------------------------
*/
   nent = Chv_countEntries(frontJ, npivot, pivotsizes, CHV_DIAGONAL) ;
   nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_BLOCK_DIAGONAL_SYM,
                                nD, nD, nent) ;
   mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
   if ( FRONTMTX_IS_SYMMETRIC(frontmtx) ) {
      SubMtx_init(mtx, frontJ->type, SUBMTX_BLOCK_DIAGONAL_SYM, 
                  J, J, nD, nD, nent) ;
   } else if ( FRONTMTX_IS_HERMITIAN(frontmtx) ) {
      SubMtx_init(mtx, frontJ->type, SUBMTX_BLOCK_DIAGONAL_HERM, 
                  J, J, nD, nD, nent) ;
   }
   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivots, &entries) ;
   IVzero(nD, pivots) ;
   IVcopy(npivot, pivots, pivotsizes) ;
   Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, entries, 
                           CHV_DIAGONAL, CHV_BY_ROWS) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n block diagonal matrix") ;
      SubMtx_writeForHumanEye(mtx, msgFile) ;
      fflush(msgFile) ;
   }
   frontmtx->nentD += nent ;
} else {
/*
   ---------------
   diagonal matrix
   ---------------
*/
   mtx = frontmtx->p_mtxDJJ[J] ;
   if ( mtx == NULL ) {
      nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_DIAGONAL, 
                                   nD, nD, nD) ;
      mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
      SubMtx_init(mtx, frontJ->type, SUBMTX_DIAGONAL, J, J, nD, nD, nD);
   }
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, entries, 
                           CHV_DIAGONAL, CHV_BY_ROWS) ;
   frontmtx->nentD += nD ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n entries in SubMtx") ;
      DVfprintf(msgFile, nent, entries) ;
      fprintf(msgFile, "\n\n block diagonal matrix") ;
      SubMtx_writeForHumanEye(mtx, msgFile) ;
      fflush(msgFile) ;
   }
}
frontmtx->p_mtxDJJ[J] = mtx ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
IVcopy(ncol, colind, colindJ) ;
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
IVcopy(nrow, rowind, rowindJ) ;
/*
   ------------------------------
   store the upper U_{J,J} matrix
   ------------------------------
*/
if ( FRONTMTX_IS_DENSE_FRONTS(frontmtx) ) {
/*
   -----------
   dense front
   -----------
*/
   nent = Chv_countEntries(frontJ, npivot, pivotsizes, 
                           CHV_STRICT_UPPER_11) ;
   if ( nent > 0 ) {
      mtx = frontmtx->p_mtxUJJ[J] ;
      if ( mtx == NULL ) {
         nbytes = SubMtx_nbytesNeeded(frontJ->type, 
                                  SUBMTX_DENSE_SUBCOLUMNS, nD, nD,nent);
         mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
         SubMtx_init(mtx, frontJ->type, SUBMTX_DENSE_SUBCOLUMNS, 
                     J, J, nD, nD, nent) ;
      }
      SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                               &firstlocs, &sizes, &entries) ;
      if ( pivotsizes != NULL ) {
         for ( jcol = ipivot = 0 ; jcol < nD ; ipivot++ ) {
            m = pivotsizes[ipivot] ;
            if ( m == 1 ) {
               firstlocs[jcol] = 0 ;
               sizes[jcol]     = jcol ;
               jcol++ ;
            } else if ( m == 2 ) {
               firstlocs[jcol] = firstlocs[jcol+1] = 0 ;
               sizes[jcol]     = sizes[jcol+1]     = jcol ;
               jcol += 2 ;
            }
         }
      } else {
         for ( jcol = 0 ; jcol < nD ; jcol++ ) {
            firstlocs[jcol] = 0 ;
            sizes[jcol]     = jcol ;
         }
      }
      Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, entries,
                              CHV_STRICT_UPPER_11, CHV_BY_COLUMNS) ;
      frontmtx->p_mtxUJJ[J] = mtx ;
      SubMtx_columnIndices(mtx, &ncol, &colind) ;
      IVcopy(ncol, colind, colindJ) ;
      SubMtx_rowIndices(mtx, &nrow, &rowind) ;
      IVcopy(nrow, rowind, rowindJ) ;
      if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
         frontmtx->nentU += nent ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n UJJ matrix") ;
         SubMtx_writeForHumanEye(mtx, msgFile) ;
         fflush(msgFile) ;
      }
   }
} else {
/*
   ------------
   sparse front
   ------------
*/
   nent = Chv_countBigEntries(frontJ, npivot, pivotsizes, 
                              CHV_STRICT_UPPER_11, droptol) ;
   if ( nent > 0 ) {
      nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_SPARSE_COLUMNS, 
                                   nD, nD, nent) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
                 "\n U_{%d,%d}, nD %d, nent %d, nbytes %d", 
                 J, J, nD, nent, nbytes) ;
         fflush(msgFile) ;
      }
      mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
      SubMtx_init(mtx, frontJ->type, SUBMTX_SPARSE_COLUMNS, 
                  J, J, nD, nD, nent) ;
      SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                             &sizes, &indices, &entries) ;
      Chv_copyBigEntriesToVector(frontJ, npivot, pivotsizes, sizes,
                                 indices, entries, CHV_STRICT_UPPER_11, 
                                 CHV_BY_COLUMNS, droptol) ;
      frontmtx->p_mtxUJJ[J] = mtx ;
      SubMtx_columnIndices(mtx, &ncol, &colind) ;
      IVcopy(ncol, colind, colindJ) ;
      SubMtx_rowIndices(mtx, &nrow, &rowind) ;
      IVcopy(nrow, rowind, rowindJ) ;
      frontmtx->nentU += nent ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n UJJ matrix") ;
         SubMtx_writeForHumanEye(mtx, msgFile) ;
         fflush(msgFile) ;
      }
   }
}
/*
   -----------------------------------
   store the upper U_{J,bnd{J}} matrix
   -----------------------------------
*/
if ( FRONTMTX_IS_DENSE_FRONTS(frontmtx) ) {
/*
   -----------
   dense front
   -----------
*/
   nent = Chv_countEntries(frontJ, npivot, pivotsizes, CHV_UPPER_12) ;
   if ( nent > 0 ) {
      mtx = frontmtx->p_mtxUJN[J] ;
      if ( mtx == NULL ) {
         nbytes = SubMtx_nbytesNeeded(frontJ->type, 
                                   SUBMTX_DENSE_COLUMNS, nD, nU, nent) ;
         mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
         SubMtx_init(mtx, frontJ->type, SUBMTX_DENSE_COLUMNS, 
                     J, nfront, nD, nU, nent) ;
      }
      SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
      Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, entries,
                              CHV_UPPER_12, CHV_BY_COLUMNS) ;
      frontmtx->p_mtxUJN[J] = mtx ;
      SubMtx_columnIndices(mtx, &ncol, &colind) ;
      IVcopy(ncol, colind, colindJ + nD) ;
      SubMtx_rowIndices(mtx, &nrow, &rowind) ;
      IVcopy(nrow, rowind, rowindJ) ;
      if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
         frontmtx->nentU += nent ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n UJJN matrix") ;
         SubMtx_writeForHumanEye(mtx, msgFile) ;
         fflush(msgFile) ;
      }
   }
} else {
/*
   ------------
   sparse front
   ------------
*/
   nent = Chv_countBigEntries(frontJ, npivot, pivotsizes, 
                              CHV_UPPER_12, droptol) ;
   if ( nent > 0 ) {
      nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_SPARSE_COLUMNS, 
                                   nD, nU, nent) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
                 "\n U_{%d,%d}, nD %d, nU %d, nent %d, nbytes %d", 
                 J, nfront, nD, nU, nent, nbytes) ;
         fflush(msgFile) ;
      }
      mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
      SubMtx_init(mtx, frontJ->type, SUBMTX_SPARSE_COLUMNS, 
               J, nfront, nD, nU, nent) ;
      SubMtx_sparseColumnsInfo(mtx, 
                             &nrow, &nent, &sizes, &indices, &entries) ;
      Chv_copyBigEntriesToVector(frontJ, npivot, pivotsizes, sizes,
                                 indices, entries, CHV_UPPER_12, 
                                 CHV_BY_COLUMNS, droptol) ;
      frontmtx->p_mtxUJN[J] = mtx ;
      SubMtx_columnIndices(mtx, &ncol, &colind) ;
      IVcopy(ncol, colind, colindJ + nD) ;
      SubMtx_rowIndices(mtx, &nrow, &rowind) ;
      IVcopy(nrow, rowind, rowindJ) ;
      frontmtx->nentU += nent ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n UJJN matrix") ;
         SubMtx_writeForHumanEye(mtx, msgFile) ;
         fflush(msgFile) ;
      }
   }
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
   ------------------------------
   store the lower L_{J,J} matrix
   ------------------------------
*/
   if ( FRONTMTX_IS_DENSE_FRONTS(frontmtx) ) {
/*
      -----------
      dense front
      -----------
*/
      nent = Chv_countEntries(frontJ, npivot, pivotsizes, 
                              CHV_STRICT_LOWER_11) ;
      if ( nent > 0 ) {
         mtx = frontmtx->p_mtxLJJ[J] ;
         if ( mtx == NULL ) {
            nbytes = SubMtx_nbytesNeeded(frontJ->type, 
                         SUBMTX_DENSE_SUBROWS, nD, nD,nent);
            mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
            SubMtx_init(mtx, frontJ->type, SUBMTX_DENSE_SUBROWS, 
                        J, J, nD, nD, nent) ;
         }
         SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                                 &firstlocs, &sizes, &entries) ;
         for ( irow = 0 ; irow < nD ; irow++ ) {
            firstlocs[irow] = 0 ;
            sizes[irow]     = irow ;
         }
         Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, 
                            entries, CHV_STRICT_LOWER_11, CHV_BY_ROWS) ;
         frontmtx->p_mtxLJJ[J] = mtx ;
         SubMtx_columnIndices(mtx, &ncol, &colind) ;
         IVcopy(ncol, colind, colindJ) ;
         SubMtx_rowIndices(mtx, &nrow, &rowind) ;
         IVcopy(nrow, rowind, rowindJ) ;
         if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
            frontmtx->nentL += nent ;
         }
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n LJJ matrix") ;
            SubMtx_writeForHumanEye(mtx, msgFile) ;
            fflush(msgFile) ;
         }
      }
   } else {
/*
      ------------
      sparse front
      ------------
*/
      nent = Chv_countBigEntries(frontJ, npivot, pivotsizes, 
                                 CHV_STRICT_LOWER_11, droptol);
      if ( nent > 0 ) {
         nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_SPARSE_ROWS, 
                                      nD, nD, nent) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    "\n L_{%d,%d}, nD %d, nent %d, nbytes %d", 
                    J, J, nD, nent, nbytes) ;
            fflush(msgFile) ;
         }
         mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
         SubMtx_init(mtx, frontJ->type, SUBMTX_SPARSE_ROWS, 
                     J, J, nD, nD, nent) ;
         SubMtx_sparseRowsInfo(mtx, 
                              &nrow, &nent, &sizes, &indices, &entries);
         Chv_copyBigEntriesToVector(frontJ, npivot, pivotsizes, sizes, 
                                  indices, entries, CHV_STRICT_LOWER_11,
                                  CHV_BY_ROWS, droptol) ;
         frontmtx->p_mtxLJJ[J] = mtx ;
         SubMtx_columnIndices(mtx, &ncol, &colind) ;
         IVcopy(ncol, colind, colindJ) ;
         SubMtx_rowIndices(mtx, &nrow, &rowind) ;
         IVcopy(nrow, rowind, rowindJ) ;
         frontmtx->nentL += nent ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n LJJ matrix") ;
            SubMtx_writeForHumanEye(mtx, msgFile) ;
            fflush(msgFile) ;
         }
      }
   }
/*
   -----------------------------------
   store the lower L_{bnd{J},J} matrix
   -----------------------------------
*/
   if ( FRONTMTX_IS_DENSE_FRONTS(frontmtx) ) {
/*
      -----------
      dense front
      -----------
*/
      nent = Chv_countEntries(frontJ, npivot, pivotsizes, CHV_LOWER_21);
      if ( nent > 0 ) {
         mtx = frontmtx->p_mtxLNJ[J] ;
         if ( mtx == NULL ) {
            nbytes = SubMtx_nbytesNeeded(frontJ->type, 
                                      SUBMTX_DENSE_ROWS, nL, nD, nent) ;
            mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
            SubMtx_init(mtx, frontJ->type, SUBMTX_DENSE_ROWS, 
                        nfront, J, nL, nD, nent) ;
         }
         SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
         Chv_copyEntriesToVector(frontJ, npivot, pivotsizes, nent, 
                                 entries, CHV_LOWER_21, CHV_BY_ROWS) ;
         frontmtx->p_mtxLNJ[J] = mtx ;
         SubMtx_columnIndices(mtx, &ncol, &colind) ;
         IVcopy(ncol, colind, colindJ) ;
         SubMtx_rowIndices(mtx, &nrow, &rowind) ;
         IVcopy(nrow, rowind, rowindJ + nD) ;
         if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
            frontmtx->nentL += nent ;
         }
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n LNJ matrix") ;
            SubMtx_writeForHumanEye(mtx, msgFile) ;
            fflush(msgFile) ;
         }
      }
   } else {
/*
      ------------
      sparse front
      ------------
*/
      nent = Chv_countBigEntries(frontJ, npivot, pivotsizes, 
                                 CHV_LOWER_21, droptol) ;
      if ( nent > 0 ) {
         nbytes = SubMtx_nbytesNeeded(frontJ->type, SUBMTX_SPARSE_ROWS, 
                                      nL, nD, nent) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, 
                    "\n L_{%d,%d}, nL %d, nD %d, nent %d, nbytes %d", 
                    nfront, J, nL, nD, nent, nbytes) ;
            fflush(msgFile) ;
         }
         mtx = SubMtxManager_newObjectOfSizeNbytes(manager, nbytes) ;
         SubMtx_init(mtx, frontJ->type, SUBMTX_SPARSE_ROWS, 
                     nfront, J, nL, nD, nent) ;
         SubMtx_sparseRowsInfo(mtx, 
                             &nrow, &nent, &sizes, &indices, &entries) ;
         Chv_copyBigEntriesToVector(frontJ, npivot, pivotsizes, sizes,
                                    indices, entries, CHV_LOWER_21, 
                                    CHV_BY_ROWS, droptol) ;
         frontmtx->p_mtxLNJ[J] = mtx ;
         SubMtx_columnIndices(mtx, &ncol, &colind) ;
         IVcopy(ncol, colind, colindJ) ;
         SubMtx_rowIndices(mtx, &nrow, &rowind) ;
         IVcopy(nrow, rowind, rowindJ + nD) ;
         frontmtx->nentL += nent ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n LNJ matrix") ;
            SubMtx_writeForHumanEye(mtx, msgFile) ;
            fflush(msgFile) ;
         }
      }
   }
}
return ; }
   
/*--------------------------------------------------------------------*/
