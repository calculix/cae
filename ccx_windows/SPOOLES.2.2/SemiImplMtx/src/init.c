/*  init.c  */

#include "../SemiImplMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- to initialize the semi-implicit matrix using as input a
              FrontMtx and a map from fronts to domains (map[J] != 0)
              or the schur complement (map[J] = 0)

   return value --
      1 -- normal return
     -1 -- semimtx is NULL
     -2 -- frontmtx is NULL
     -3 -- inpmtx is NULL
     -4 -- frontmapIV is NULL
     -5 -- frontmapIV is invalid
     -6 -- unable to create domains' front matrix
     -7 -- unable to create schur complement front matrix

   created -- 98oct17, cca
   ------------------------------------------------------------------
*/
int
SemiImplMtx_initFromFrontMtx (
   SemiImplMtx   *semimtx,
   FrontMtx      *frontmtx,
   InpMtx        *inpmtx,
   IV            *frontmapIV,
   int           msglvl,
   FILE          *msgFile
) {
FrontMtx   *domMtx, *schurMtx ;
InpMtx     *A12, *A21 ;
int        ii, J, ncol, nfront, nrow, rc, size ;
int        *cols, *frontmap, *rows ;
IV         *domColsIV, *domidsIV, *domRowsIV, 
           *schurColsIV, *schuridsIV, *schurRowsIV ;
/*
   --------------
   check the data
   --------------
*/
if ( semimtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n semimtx is NULL\n") ;
   return(-1) ;
}
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n frontmtx is NULL\n") ;
   return(-2) ;
}
if ( inpmtx == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n inpmtx is NULL\n") ;
   return(-3) ;
}
if ( frontmapIV == NULL ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n frontmapIV is NULL\n") ;
   return(-4) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
IV_sizeAndEntries(frontmapIV, &size, &frontmap) ;
if ( nfront != size ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n nfront %d, size of front map %d\n", nfront, size) ;
   return(-5) ;
}
domidsIV   = IV_new() ;
schuridsIV = IV_new() ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( frontmap[J] == 0 ) {
      IV_push(schuridsIV, J) ;
   } else if ( frontmap[J] > 0 ) {
      IV_push(domidsIV, J) ;
   } else {
      fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
              "\n frontmap[%d] = %d, invalid\n", J, frontmap[J]) ;
      IV_free(domidsIV) ;
      IV_free(schuridsIV) ;
      return(-5) ;
   }
}
/*
   -----------------------------------------------------------
   clear the data for the semi-implicit matrix and set scalars
   -----------------------------------------------------------
*/
SemiImplMtx_clearData(semimtx) ;
semimtx->neqns = frontmtx->neqns ;
semimtx->type  = frontmtx->type  ;
semimtx->symmetryflag = frontmtx->symmetryflag ;
/*
   ----------------------------------------------
   get the front matrix that contains the domains
   ----------------------------------------------
*/
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n working on domain front matrix") ;
   fflush(msgFile) ;
}
domMtx = semimtx->domainMtx = FrontMtx_new() ;
domRowsIV = semimtx->domRowsIV = IV_new() ;
domColsIV = semimtx->domColsIV = IV_new() ;
rc = FrontMtx_initFromSubmatrix(domMtx, frontmtx, domidsIV, 
                                domRowsIV, domColsIV, msglvl, msgFile) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n unable to initialize the domains' front matrix"
           "\n error return = %d\n", rc) ;
   return(-6) ;
}
semimtx->ndomeqns = IV_size(domRowsIV) ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n---------------------------------------- ") ;
   fprintf(msgFile, "\n\n submatrix for domains") ;
   FrontMtx_writeForHumanEye(domMtx, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   FrontMtx_writeForMatlab(domMtx, "L11", "D11", "U11", msgFile) ;
   IV_writeForMatlab(domRowsIV, "domrows", msgFile) ;
   IV_writeForMatlab(domColsIV, "domcols", msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------
   get the front matrix that contains the schur complement
   -------------------------------------------------------
*/
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n working on domain front matrix") ;
   fflush(msgFile) ;
}
schurMtx = semimtx->schurMtx = FrontMtx_new() ;
schurRowsIV = semimtx->schurRowsIV = IV_new() ;
schurColsIV = semimtx->schurColsIV = IV_new() ;
rc = FrontMtx_initFromSubmatrix(schurMtx, frontmtx, schuridsIV, 
                            schurRowsIV, schurColsIV, msglvl, msgFile) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n unable to initialize the schur complement front matrix"
           "\n error return = %d\n", rc) ;
   return(-6) ;
}
semimtx->nschureqns = IV_size(schurRowsIV) ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n---------------------------------------- ") ;
   fprintf(msgFile, "\n\n submatrix for schur complement") ;
   FrontMtx_writeForHumanEye(schurMtx, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   FrontMtx_writeForMatlab(schurMtx, "L22", "D22", "U22", msgFile) ;
   IV_writeForMatlab(schurRowsIV, "schurrows", msgFile) ;
   IV_writeForMatlab(schurColsIV, "schurcols", msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   get the A12 InpMtx object
   -------------------------
*/
A12 = semimtx->A12 = InpMtx_new() ;
rc = InpMtx_initFromSubmatrix(A12, inpmtx, domRowsIV, schurColsIV,
                              semimtx->symmetryflag, msglvl, msgFile) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
           "\n unable to create A21 matrix"
           "\n error return = %d\n", rc) ;
   return(-6) ;
}
InpMtx_changeCoordType(A12, INPMTX_BY_ROWS) ;
InpMtx_changeStorageMode(A12, INPMTX_BY_VECTORS) ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n---------------------------------------- ") ;
   fprintf(msgFile, "\n\n domRowsIV ") ;
   IV_writeForHumanEye(domRowsIV, msgFile) ;
   fprintf(msgFile, "\n\n schurColsIV ") ;
   IV_writeForHumanEye(schurColsIV, msgFile) ;
   fprintf(msgFile, "\n\n A12 matrix") ;
   InpMtx_writeForHumanEye(A12, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n A12 = zeros(%d,%d) ;",
           IV_size(domRowsIV), IV_size(schurColsIV)) ;
   InpMtx_writeForMatlab(A12, "A12", msgFile) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
   -------------------------
   get the A21 InpMtx object
   -------------------------
*/
   A21 = semimtx->A21 = InpMtx_new() ;
   rc = InpMtx_initFromSubmatrix(A21, inpmtx, schurRowsIV, domColsIV,
                              semimtx->symmetryflag, msglvl, msgFile) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n error in SemiImplMtx_initFromFrontMtx()"
              "\n unable to create A21 matrix"
              "\n error return = %d\n", rc) ;
      return(-6) ;
   }
   InpMtx_changeCoordType(A21, INPMTX_BY_COLUMNS) ;
   InpMtx_changeStorageMode(A21, INPMTX_BY_VECTORS) ;
   if ( msglvl > 4 ) {
      fprintf(msgFile, "\n\n--------------------------------------- ") ;
      fprintf(msgFile, "\n\n schurRowsIV ") ;
      IV_writeForHumanEye(schurRowsIV, msgFile) ;
      fprintf(msgFile, "\n\n domColsIV ") ;
      IV_writeForHumanEye(domColsIV, msgFile) ;
      fprintf(msgFile, "\n\n A21 matrix") ;
      InpMtx_writeForHumanEye(A21, msgFile) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 4 ) {
      fprintf(msgFile, "\n\n A21 = zeros(%d,%d) ;",
              IV_size(schurRowsIV), IV_size(domColsIV)) ;
      InpMtx_writeForMatlab(A21, "A21", msgFile) ;
      fflush(msgFile) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(domidsIV) ;
IV_free(schuridsIV) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   purpose -- to fill submtx with a submatrix of the front matrix.
      the fronts that form the submatrix are found in frontidsIV.

      all information in submtx is local, front #'s are from 0 to
      one less than the number of fronts in the submatrix, equation
      #'s are from 0 to one less than the number of rows and columns
      in the submatrix. the global row and column ids for the submatrix
      are stored in rowsIV and colsIV on return.

   return values ---
      1 -- normal return
     -1 -- submtx is NULL
     -2 -- frontmtx is NULL
     -3 -- frontmtx is not in 2-D mode
     -4 -- frontidsIV is NULL
     -5 -- frontidsIV is invalid
     -6 -- rowsIV is NULL
     -7 -- colsIV is NULL
     -8 -- unable to create front tree
     -9 -- unable to create symbfacIVL
    -10 -- unable to create coladjIVL
    -11 -- unable to create rowadjIVL
    -12 -- unable to create upperblockIVL
    -13 -- unable to create lowerblockIVL

   created -- 98oct17, cca
   --------------------------------------------------------------------
*/
int
FrontMtx_initFromSubmatrix (
   FrontMtx   *submtx,
   FrontMtx   *frontmtx,
   IV         *frontidsIV,
   IV         *rowsIV,
   IV         *colsIV,
   int        msglvl,
   FILE       *msgFile
) {
ETree    *etreeSub ;
int      ii, J, Jsub, K, Ksub, ncol, nfront, nfrontSub, neqnSub, nJ,
         nrow, offset, rc, size, vSub ;
int      *bndwghts, *colind, *colmap, *cols, *frontSubIds, 
         *list, *nodwghts, *rowind, *rowmap, *rows ;
IV       *frontsizesIVsub, *vtxIV ;
IVL      *coladjIVLsub, *lowerblockIVLsub, *rowadjIVLsub, 
         *symbfacIVLsub, *upperblockIVLsub ;
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( submtx == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n submtx is NULL\n") ;
   return(-1) ;
}
if ( frontmtx == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n frontmtx is NULL\n") ;
   return(-2) ;
}
if ( ! FRONTMTX_IS_2D_MODE(frontmtx) ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n frontmtx mode is not 2D\n") ;
   return(-3) ;
}
if ( frontidsIV == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n frontidsIV is NULL\n") ;
   return(-4) ;
}
nfront = FrontMtx_nfront(frontmtx) ;
IV_sizeAndEntries(frontidsIV, &nfrontSub, &frontSubIds) ;
if ( nfrontSub < 0 || nfrontSub > nfront ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n invalid frontidsIV"
           "\n nfrontSub = %d, nfront %d\n", nfrontSub, nfront) ;
   return(-5) ;
}
for ( ii = 0 ; ii < nfrontSub ; ii++ ) {
   if ( (J = frontSubIds[ii]) < 0 || J >= nfront ) {
      fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
              "\n invalid frontidsIV"
              "\n frontSubIds[%d] = %d, nfront = %d\n",
              ii, J, nfront) ;
      return(-5) ;
   }
}
if ( rowsIV == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n rowsIV is NULL\n") ;
   return(-6) ;
}
if ( colsIV == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n colsIV is NULL\n") ;
   return(-7) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   clear the data for the submatrix and set the 
   scalar values (some inherited from the global matrix)
   -----------------------------------------------------
*/
FrontMtx_clearData(submtx) ;
submtx->nfront       = nfrontSub ;
submtx->type         = frontmtx->type ;
submtx->symmetryflag = frontmtx->symmetryflag ;
submtx->sparsityflag = frontmtx->sparsityflag ;
submtx->pivotingflag = frontmtx->pivotingflag ;
submtx->dataMode     = FRONTMTX_2D_MODE ;
/*
   ---------------------------------------------------------------
   initialize the front tree for the submatrix.

   note: on return, vtxIV is filled with the vertices originally
   in the submatrix, (pivoting may change this), needed to find
   symbolic factorization IVL object

   note: at return, the boundary weights are likely to be invalid,
   since we have no way of knowing what boundary indices for a
   front are really in the domain. this will be changed after we
   have the symbolic factorization.
   ---------------------------------------------------------------
*/
etreeSub = submtx->frontETree = ETree_new() ;
vtxIV = IV_new() ;
rc = ETree_initFromSubtree(etreeSub, frontidsIV, 
                           frontmtx->frontETree, vtxIV) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
         "\n unable to create submatrix's front ETree, rc = %d\n", rc) ;
   return(-8) ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n submatrix ETree") ;
   ETree_writeForHumanEye(etreeSub, msgFile) ;
   fprintf(msgFile, "\n\n submatrix original equations") ;
   IV_writeForHumanEye(vtxIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------------------
   set the # of equations (perhap temporarily if pivoting 
   has delayed some rows and columns), and the tree.
   ------------------------------------------------------
*/
submtx->neqns = neqnSub = IV_size(vtxIV) ;
submtx->tree  = etreeSub->tree ;
/*
   -----------------------------------------------------
   initialize the symbolic factorization for the subtree
   -----------------------------------------------------
*/
symbfacIVLsub = submtx->symbfacIVL = IVL_new() ;
rc = IVL_initFromSubIVL(symbfacIVLsub, frontmtx->symbfacIVL,
                        frontidsIV, vtxIV) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
         "\n unable to create submatrix's symbfac, rc = %d\n", rc) ;
   return(-9) ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n submatrix symbolic factorizatio") ;
   IVL_writeForHumanEye(symbfacIVLsub, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   adjust the boundary weights of the front tree
   ---------------------------------------------
*/
nodwghts = ETree_nodwghts(etreeSub) ;
bndwghts = ETree_bndwghts(etreeSub) ;
for ( J = 0 ; J < nfrontSub ; J++ ) {
   IVL_listAndSize(symbfacIVLsub, J, &size, &list) ;
   bndwghts[J] = size - nodwghts[J] ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n submatrix ETree after bndweight adjustment") ;
   ETree_writeForHumanEye(etreeSub, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------
   set the front sizes for the submatrix
   -------------------------------------
*/
frontsizesIVsub = submtx->frontsizesIV = IV_new() ;
IV_init(frontsizesIVsub, nfrontSub, NULL) ;
IVgather(nfrontSub, IV_entries(frontsizesIVsub), 
         IV_entries(frontmtx->frontsizesIV),
         IV_entries(frontidsIV)) ;
neqnSub = submtx->neqns = IV_sum(frontsizesIVsub) ;
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n %d equations in submatrix", neqnSub) ;
   fprintf(msgFile, "\n\n front sizes for submatrix") ;
   IV_writeForHumanEye(frontsizesIVsub, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------------------
   fill rowsIV and colsIV with the row and column ids of the submatrix
   -------------------------------------------------------------------
*/
IV_setSize(rowsIV, neqnSub) ;
IV_setSize(colsIV, neqnSub) ;
rows = IV_entries(rowsIV) ;
cols = IV_entries(colsIV) ;
for ( Jsub = offset = 0 ; Jsub < nfrontSub ; Jsub++ ) {
   if ( (nJ = FrontMtx_frontSize(submtx, Jsub)) > 0 ) {
      J = frontSubIds[Jsub] ;
      FrontMtx_columnIndices(frontmtx, J, &size, &list) ;
      IVcopy(nJ, cols + offset, list) ;
      FrontMtx_rowIndices(frontmtx, J, &size, &list) ;
      IVcopy(nJ, rows + offset, list) ;
      offset += nJ ;
   }
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n row ids for submatrix") ;
   IV_writeForHumanEye(rowsIV, msgFile) ;
   fprintf(msgFile, "\n\n column ids for submatrix") ;
   IV_writeForHumanEye(colsIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   get the row and column adjacencies
   ----------------------------------
*/
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   submtx->neqns = neqnSub ;
   coladjIVLsub  = submtx->coladjIVL = IVL_new() ;
   rc = IVL_initFromSubIVL(coladjIVLsub, frontmtx->coladjIVL,
                           frontidsIV, colsIV) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n unable to create submatrix's coladjIVL, rc = %d\n", rc) ;
      return(-10) ;
   }
   if ( msglvl > 4 ) {
      fprintf(msgFile, "\n\n submatrix col adjacency") ;
      IVL_writeForHumanEye(coladjIVLsub, msgFile) ;
      fflush(msgFile) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      rowadjIVLsub = submtx->rowadjIVL = IVL_new() ;
      rc = IVL_initFromSubIVL(rowadjIVLsub, frontmtx->rowadjIVL,
                              frontidsIV, rowsIV) ;
      if ( rc != 1 ) {
         fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n unable to create submatrix's rowadjIVL, rc = %d\n", rc) ;
         return(-11) ;
      }
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n\n submatrix row adjacency") ;
         IVL_writeForHumanEye(rowadjIVLsub, msgFile) ;
         fflush(msgFile) ;
      }
   }
}
IV_free(vtxIV) ;
/*
   ----------------------------------------------
   get the rowmap[] and colmap[] vectors,
   needed to translate indices in the submatrices
   ----------------------------------------------
*/
colmap = IVinit(frontmtx->neqns, -1) ;
for ( ii = 0 ; ii < neqnSub ; ii++ ) {
   colmap[cols[ii]] = ii ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   rowmap = IVinit(frontmtx->neqns, -1) ;
   for ( ii = 0 ; ii < neqnSub ; ii++ ) {
      rowmap[rows[ii]] = ii ;
   }
} else {
   rowmap = colmap ;
}
/*
   -----------------------------------------------------------
   get the upper and lower block IVL objects for the submatrix
   -----------------------------------------------------------
*/
upperblockIVLsub = submtx->upperblockIVL = IVL_new() ;
rc = IVL_initFromSubIVL(upperblockIVLsub, frontmtx->upperblockIVL,
                        frontidsIV, frontidsIV) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
        "\n unable to create upperblockIVL, rc = %d\n", rc) ;
   return(-12) ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n\n upper block adjacency IVL object") ;
   IVL_writeForHumanEye(upperblockIVLsub, msgFile) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   lowerblockIVLsub = submtx->lowerblockIVL = IVL_new() ;
   rc = IVL_initFromSubIVL(lowerblockIVLsub, frontmtx->lowerblockIVL,
                           frontidsIV, frontidsIV) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_initFromSubmatrix()"
           "\n unable to create lowerblockIVL, rc = %d\n", rc) ;
      return(-13) ;
   }
   if ( msglvl > 4 ) {
      fprintf(msgFile, "\n\n lower block adjacency IVL object") ;
      IVL_writeForHumanEye(lowerblockIVLsub, msgFile) ;
      fflush(msgFile) ;
   }
}
/*
   ----------------------------------------------------------------
   allocate the vector and hash table(s) for the factor submatrices
   ----------------------------------------------------------------
*/
ALLOCATE(submtx->p_mtxDJJ, struct _SubMtx *, nfrontSub) ;
for ( J = 0 ; J < nfrontSub ; J++ ) {
   submtx->p_mtxDJJ[J] = NULL ;
}
submtx->upperhash = I2Ohash_new() ;
I2Ohash_init(submtx->upperhash, nfrontSub, nfrontSub, nfrontSub) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   submtx->lowerhash = I2Ohash_new() ;
   I2Ohash_init(submtx->lowerhash, nfrontSub, nfrontSub, nfrontSub) ;
}
/*
   -----------------------------------------------------------------
   remove the diagonal submatrices from the factor matrix
   and insert into the submatrix object. note: front row and column
   ids must be changed to their local values, and the row and column
   indices must be mapped to local indices.
   -----------------------------------------------------------------
*/
for ( Jsub = 0 ; Jsub < nfrontSub ; Jsub++ ) {
   J = frontSubIds[Jsub] ;
   if ( (mtx = frontmtx->p_mtxDJJ[J]) != NULL ) {
      SubMtx_setIds(mtx, Jsub, Jsub) ;
      SubMtx_columnIndices(mtx, &ncol, &colind) ;
      IVgather(ncol, colind, colmap, colind) ;
      SubMtx_rowIndices(mtx, &nrow, &rowind) ;
      IVgather(nrow, rowind, rowmap, rowind) ;
      submtx->p_mtxDJJ[Jsub] = mtx ;
      frontmtx->p_mtxDJJ[J]  = NULL ;
      submtx->nentD += mtx->nent ;
   }
}
/*
   ----------------------------------------------------------------
   remove the upper triangular submatrices from the factor matrix
   and insert into the submatrix object. note: front row and column
   ids must be changed to their local values. if the matrix is on
   the diagonal, i.e., U(J,J), its row and column indices must be 
   mapped to local indices.
   ----------------------------------------------------------------
*/
for ( Jsub = 0 ; Jsub < nfrontSub ; Jsub++ ) {
   J = frontSubIds[Jsub] ;
   FrontMtx_upperAdjFronts(submtx, Jsub, &size, &list) ;
   for ( ii = 0 ; ii < size ; ii++ ) {
      Ksub = list[ii] ;
      K = frontSubIds[Ksub] ;
      if ( 1 == I2Ohash_remove(frontmtx->upperhash, 
                               J, K, (void *) &mtx) ) {
         SubMtx_setIds(mtx, Jsub, Ksub) ;
         if ( K == J ) {
            SubMtx_columnIndices(mtx, &ncol, &colind) ;
            IVgather(ncol, colind, colmap, colind) ;
            SubMtx_rowIndices(mtx, &nrow, &rowind) ;
            IVgather(nrow, rowind, rowmap, rowind) ;
         }
         I2Ohash_insert(submtx->upperhash, Jsub, Ksub, (void *) mtx) ;
         submtx->nentU += mtx->nent ;
      }
   }
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
   ----------------------------------------------------------------
   remove the lower triangular submatrices from the factor matrix
   and insert into the submatrix object. note: front row and column
   ids must be changed to their local values. if the matrix is on
   the diagonal, i.e., L(J,J), its row and column indices must be 
   mapped to local indices.
   ----------------------------------------------------------------
*/
   for ( Jsub = 0 ; Jsub < nfrontSub ; Jsub++ ) {
      J = frontSubIds[Jsub] ;
      FrontMtx_lowerAdjFronts(submtx, Jsub, &size, &list) ;
      for ( ii = 0 ; ii < size ; ii++ ) {
         Ksub = list[ii] ;
         K = frontSubIds[Ksub] ;
         if ( 1 == I2Ohash_remove(frontmtx->lowerhash, 
                                  K, J, (void *) &mtx) ) {
            SubMtx_setIds(mtx, Ksub, Jsub) ;
            if ( K == J ) {
               SubMtx_columnIndices(mtx, &ncol, &colind) ;
               IVgather(ncol, colind, colmap, colind) ;
               SubMtx_rowIndices(mtx, &nrow, &rowind) ;
               IVgather(nrow, rowind, rowmap, rowind) ;
            }
            I2Ohash_insert(submtx->lowerhash, Ksub, Jsub, (void *) mtx);
            submtx->nentL += mtx->nent ;
         }
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(colmap) ;
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   IVfree(rowmap) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
