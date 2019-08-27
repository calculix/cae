/*  instance.c  */

#include "../SubMtx.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   purpose -- fill *prowid with the row id
              and  *pcolid with the column id

   created -- 98may01, cca
   ------------------------------------------
*/
void
SubMtx_ids (
   SubMtx   *mtx,
   int      *prowid,
   int      *pcolid
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || prowid == NULL || pcolid == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_ids(%p,%p,%p)"
           "\n bad input\n", mtx, prowid, pcolid) ;
   exit(-1) ;
}
*prowid = mtx->rowid ;
*pcolid = mtx->colid ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- set the row and column ids

   created -- 98may01, cca
   -------------------------------------
*/
void
SubMtx_setIds (
   SubMtx   *mtx,
   int      rowid,
   int      colid
) {
int   *buffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_ids(%p,%d,%d)"
           "\n bad input\n", mtx, rowid, colid) ;
   exit(-1) ;
}
buffer = (int *) mtx->wrkDV.vec ;
buffer[2] = mtx->rowid = rowid ;
buffer[3] = mtx->colid = colid ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- fill *pnrow with the # of rows
              and  *pncol with the # of columns
              and  *pnent with the # of matrix entries

   created -- 98may01, cca
   ---------------------------------------------------
*/
void
SubMtx_dimensions (
   SubMtx   *mtx,
   int      *pnrow,
   int      *pncol,
   int      *pnent
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || pnrow == NULL || pncol == NULL || pnent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_ids(%p,%p,%p,%p)"
           "\n bad input\n", mtx, pnrow, pncol, pnent) ;
   exit(-1) ;
}
*pnrow = mtx->nrow ;
*pncol = mtx->ncol ;
*pnent = mtx->nent ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- 
     fill *pnrow with the number of rows
     if prowind != NULL then
        fill *prowind with the base location of the row indices
     endif

   created -- 98may01, cca
   ------------------------------------------------------------
*/
void
SubMtx_rowIndices (
   SubMtx   *mtx,
   int      *pnrow,
   int      **prowind
) {
*pnrow = mtx->nrow ;
if ( prowind != NULL ) {
   *prowind = (int *) mtx->wrkDV.vec + 7 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- 
     fill *pncol with the number of column
     if pcolind != NULL then
        fill *prowind with the base location of the column indices
     endif

   created -- 98may01, cca
   ---------------------------------------------------------------
*/
void
SubMtx_columnIndices (
   SubMtx   *mtx,
   int      *pncol,
   int      **pcolind
) {
*pncol = mtx->ncol ;
if ( pcolind != NULL ) {
   *pcolind = (int *) mtx->wrkDV.vec + 7 + mtx->nrow ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   purpose -- for dense storage
     *pnrow    with mtx->nrow
     *pncol    with mtx->ncol
     *pinc1    with row increment
     *pinc2    with column increment
     *pentries with mtx->entries

   created -- 98may01, cca
   ---------------------------------
*/
void
SubMtx_denseInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pncol,
   int        *pinc1,
   int        *pinc2,
   double     **pentries
) {
double   *dbuffer ;
int      nint ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || pnrow == NULL || pncol == NULL 
   || pinc1 == NULL || pinc2 == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_denseInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           mtx, pnrow, pncol, pinc1, pinc2, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_denseInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           mtx, pnrow, pncol, pinc1, pinc2, pentries, mtx->type) ;
            exit(-1) ;
}
if ( ! (SUBMTX_IS_DENSE_ROWS(mtx) || SUBMTX_IS_DENSE_COLUMNS(mtx)) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_denseInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad mode %d"
           "\n must be SUBMTX_DENSE_ROWS or SUBMTX_DENSE_COLUMNS\n",
           mtx, pnrow, pncol, pinc1, pinc2, pentries, mtx->mode) ;
   exit(-1) ;
}
*pnrow = mtx->nrow ;
*pncol = mtx->ncol ;
if ( SUBMTX_IS_DENSE_ROWS(mtx) ) {
   *pinc1 = mtx->ncol ;
   *pinc2 = 1 ;
} else {
   *pinc1 = 1 ;
   *pinc2 = mtx->nrow ;
}
dbuffer = mtx->wrkDV.vec ;
nint = 7 + mtx->nrow + mtx->ncol ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   purpose -- for sparse rows, fill
     *pnrow    with # of rows
     *pnent    with # of indices
     *psizes   with sizes[] of rows
     *pindices with indices[] for column indices
     *pentries with entries[] for matrix entries

   created -- 98may01, cca
   ---------------------------------------------
*/
void
SubMtx_sparseRowsInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pnent,
   int        **psizes,
   int        **pindices,
   double     **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pnrow == NULL  || pnent == NULL 
   || psizes == NULL || pindices == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_sparseRowsInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           mtx, pnrow, pnent, psizes, pindices, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_sparseRowsInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
           mtx, pnrow, pnent, psizes, pindices, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_SPARSE_ROWS(mtx) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseRowsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad mode %d, must be SUBMTX_SPARSE_ROWS\n",
         mtx, pnrow, pnent, psizes, pindices, pentries, mtx->mode) ;
   exit(-1) ;
}
*pnrow    = mtx->nrow ;
*pnent    = mtx->nent ;
dbuffer   = mtx->wrkDV.vec ;
ibuffer   = (int *) dbuffer ;
nint      = 7 + mtx->nrow + mtx->ncol ;
*psizes   = ibuffer + nint ;
nint      += mtx->nrow ;
*pindices = ibuffer + nint ;
nint      += mtx->nent ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- for sparse columns, fill
     *pncol    with # of columns
     *pnent    with # of matrix entries
     *psizes   with sizes[ncol], column sizes
     *pindices with indices[nent], matrix row ids
     *pentries with entries[nent], matrix entries

   created -- 98may01, cca
   ----------------------------------------------
*/
void
SubMtx_sparseColumnsInfo (
   SubMtx     *mtx,
   int        *pncol,
   int        *pnent,
   int        **psizes,
   int        **pindices,
   double     **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pncol == NULL  || pnent == NULL 
   || psizes == NULL || pindices == NULL || pentries == NULL ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseColumnsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad input\n",
         mtx, pncol, pnent, psizes, pindices, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseColumnsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
         mtx, pncol, pnent, psizes, pindices, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_SPARSE_COLUMNS(mtx) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseColumnsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad mode %d"
         "\n must be SUBMTX_SPARSE_COLUMNS\n",
         mtx, pncol, pnent, psizes, pindices, pentries, mtx->mode) ;
   exit(-1) ;
}
*pncol    = mtx->ncol ;
*pnent    = mtx->nent ;
dbuffer   = mtx->wrkDV.vec ;
ibuffer   = (int *) dbuffer ;
nint      = 7 + mtx->nrow + mtx->ncol ;
*psizes   = ibuffer + nint ;
nint      += mtx->ncol ;
*pindices = ibuffer + nint ;
nint      += mtx->nent ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- for sparse triples, fill
     *pnent    with # of matrix entries
     *prowids  with rowids[nent], row ids of entries
     *prowids  with colids[nent], column ids of entries
     *pentries with entries[nent], matrix entries

   created -- 98may01, cca
   ----------------------------------------------------
*/
void
SubMtx_sparseTriplesInfo (
   SubMtx     *mtx,
   int        *pnent,
   int        **prowids,
   int        **pcolids,
   double     **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pnent == NULL 
   || prowids == NULL || pcolids == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_sparseTriplesInfo(%p,%p,%p,%p,%p)"
           "\n bad input\n",
           mtx, pnent, prowids, pcolids, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseTriplesInfo(%p,%p,%p,%p,%p)"
         "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
         mtx, pnent, prowids, pcolids, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_SPARSE_TRIPLES(mtx) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_sparseTriplesInfo(%p,%p,%p,%p,%p)"
         "\n bad mode %d"
         "\n must be SUBMTX_SPARSE_TRIPLES\n",
         mtx, pnent, prowids, pcolids, pentries, mtx->mode) ;
   exit(-1) ;
}
*pnent   = mtx->nent ;
dbuffer  = mtx->wrkDV.vec ;
ibuffer  = (int *) dbuffer ;
nint     = 7 + mtx->nrow + mtx->ncol ;
*prowids = ibuffer + nint ;
nint     += mtx->nent ;
*pcolids = ibuffer + nint ;
nint     += mtx->nent ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- for dense subrows, fill
     *pnrow      with # of rows
     *pnent      with # of matrix entries
     *pfirstlocs with firstlocs[nrow], column of first nonzero
     *psizes     with sizes[nrow], number of nonzero columns
     *pentries   with entries[nent], matrix entries

   created -- 98may01, cca
   -----------------------------------------------------------
*/
void
SubMtx_denseSubrowsInfo (
   SubMtx     *mtx,
   int        *pnrow,
   int        *pnent,
   int        **pfirstlocs,
   int        **psizes,
   double     **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pnrow == NULL || pnent == NULL
   || pfirstlocs == NULL || psizes == NULL || pentries == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_denseSubrowsInfo(%p,%p,%p,%p,%p,%p)"
           "\n bad input\n",
           mtx, pnrow, pnent, pfirstlocs, psizes, pentries) ;
   if ( mtx != NULL ) {
      SubMtx_writeForHumanEye(mtx, stderr) ;
   }
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_denseSubrowsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
         mtx, pnrow, pnent, pfirstlocs, psizes, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_SUBROWS(mtx) ) {
   fprintf(stderr, 
         "\n fatal error in SubMtx_denseSubrowsInfo(%p,%p,%p,%p,%p,%p)"
         "\n bad mode %d"
         "\n must be SUBMTX_DENSE_SUBROWS\n",
         mtx, pnrow, pnent, pfirstlocs, psizes, pentries, mtx->mode) ;
   exit(-1) ;
}
*pnrow  = mtx->nrow ;
*pnent  = mtx->nent ;
dbuffer = mtx->wrkDV.vec ;
ibuffer = (int *) dbuffer ;
nint    = 7 + mtx->nrow + mtx->ncol ;
*pfirstlocs = ibuffer + nint ;
nint    += mtx->nrow ;
*psizes = ibuffer + nint ;
nint    += mtx->nrow ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- for dense subcolumns, fill
     *pncol      with # of columns
     *pnent      with # of matrix entries
     *pfirstlocs with firstlocs[ncol], row of first nonzero
     *psizes     with sizes[ncol], number of nonzero rows
     *pentries   with entries[nent], matrix entries

   created -- 98may01, cca
   -----------------------------------------------------------
*/
void
SubMtx_denseSubcolumnsInfo (
   SubMtx   *mtx,
   int      *pncol,
   int      *pnent,
   int      **pfirstlocs,
   int      **psizes,
   double   **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL 
   || pfirstlocs == NULL || psizes == NULL || pentries == NULL ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_denseSubcolumnsInfo(%p,%p,%p,%p,%p,%p)"
        "\n bad input\n",
        mtx, pncol, pnent, pfirstlocs, psizes, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_denseSubcolumsInfo(%p,%p,%p,%p,%p,%p)"
        "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
        mtx, pncol, pnent, pfirstlocs, psizes, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DENSE_SUBCOLUMNS(mtx) ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_denseSubcolumnsInfo(%p,%p,%p,%p,%p,%p)"
       "\n bad mode %d"
       "\n must be SUBMTX_DENSE_SUBCOLUMNS\n",
       mtx, pncol, pnent, pfirstlocs, psizes, pentries, mtx->mode) ;
   exit(-1) ;
}
*pncol = mtx->ncol ;
*pnent = mtx->nent ;
dbuffer = mtx->wrkDV.vec ;
ibuffer = (int *) dbuffer ;
nint = 7 + mtx->nrow + mtx->ncol ;
*pfirstlocs = ibuffer + nint ;
nint += mtx->ncol ;
*psizes = ibuffer + nint ;
nint += mtx->ncol ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- for a diagonal matrix, fill
     *pncol      with # of columns
     *pentries   with entries[nent], matrix entries

   created -- 98may01, cca
   ------------------------------------------------
*/
void
SubMtx_diagonalInfo (
   SubMtx   *mtx,
   int      *pncol,
   double   **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || pncol == NULL || pentries == NULL ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_diagonalInfo(%p,%p,%p)"
        "\n bad input\n",
        mtx, pncol, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_diagonalInfo(%p,%p,%p)"
        "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
        mtx, pncol, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_DIAGONAL(mtx) ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_diagonalInfo(%p,%p,%p)"
       "\n bad mode %d"
       "\n must be SUBMTX_DIAGONAL\n",
       mtx, pncol, pentries, mtx->mode) ;
   exit(-1) ;
}
*pncol = mtx->ncol ;
dbuffer = mtx->wrkDV.vec ;
ibuffer = (int *) dbuffer ;
nint = 7 + mtx->nrow + mtx->ncol ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- for a block diagonal symmetric matrix, fill
     *pncol        with # of columns
     *pnent        with # of entries
     *ppivotsizes  with pivotsizes[ncol]
     *pentries     with entries[nent], matrix entries

   created -- 98may01, cca
   ------------------------------------------------------
*/
void
SubMtx_blockDiagonalInfo (
   SubMtx   *mtx,
   int      *pncol,
   int      *pnent,
   int      **ppivotsizes,
   double   **pentries
) {
double   *dbuffer ;
int      nint ;
int      *ibuffer ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL 
   || pncol == NULL || pnent == NULL 
   || ppivotsizes == NULL || pentries == NULL ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_blockDiagonalInfo(%p,%p,%p,%p,%p)"
        "\n bad input\n",
        mtx, pncol, pnent, ppivotsizes, pentries) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtx) || SUBMTX_IS_COMPLEX(mtx)) ) {
   fprintf(stderr, 
        "\n fatal error in SubMtx_blockDiagonalInfo(%p,%p,%p,%p,%p)"
        "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n",
        mtx, pncol, pnent, ppivotsizes, pentries, mtx->type) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_BLOCK_DIAGONAL_SYM(mtx)
        || SUBMTX_IS_BLOCK_DIAGONAL_HERM(mtx)) ) {
   fprintf(stderr, 
"\n fatal error in SubMtx_blockDiagonalInfo(%p,%p,%p,%p,%p)"
"\n bad mode %d"
"\n must be SUBMTX_BLOCK_DIAGONAL_SYM or SUBMTX_BLOCK_DIAGONAL_HERM \n",
mtx, pncol, pnent, ppivotsizes, pentries, mtx->mode) ;
   exit(-1) ;
}
*pncol = mtx->ncol ;
*pnent = mtx->nent ;
dbuffer = mtx->wrkDV.vec ;
ibuffer = (int *) dbuffer ;
nint = 7 + 2*mtx->nrow ;
*ppivotsizes = ibuffer + nint ;
nint += mtx->nrow ;
if ( sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + nint ;
} else if ( 2*sizeof(int) == sizeof(double) ) {
   *pentries = dbuffer + (nint+1)/2 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to find matrix entry (irow,jcol) if present.

   return value --
     if entry (irow,jcol) is not present then
        *pValue is 0.0
        return value is -1
     else entry (irow,jcol) is present then
        *pValue is the matrix entry
        return value is offset into entries array 
     endif

   created -- 98may01, cca
   -------------------------------------------------------
*/
int
SubMtx_realEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   *pValue
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || irow < 0 || irow >= mtx->nrow || jcol < 0 
   || jcol >= mtx->ncol || pValue == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_realEntry(%p,%d,%d,%p)"
           "\n bad input\n", mtx, irow, jcol, pValue) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_REAL(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_realEntry(%p,%d,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, irow, jcol, pValue, mtx->type) ;
   exit(-1) ;
}
*pValue = 0 ;
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, ncol, nrow, offset ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( irow < 0 || irow >= nrow || jcol < 0 || jcol >= ncol ) {
      return(-1) ;
   }
   offset  = irow*inc1 + jcol*inc2 ;
   *pValue = entries[offset] ;
   return(offset) ;
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jj, nent, nrow, offset, *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   if ( irow < 0 || irow >= nrow ) {
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < irow ; ii++ ) {
      offset += sizes[ii] ;
   }
   for ( ii = 0, jj = offset ; ii < sizes[irow] ; ii++, jj++ ) {
      if ( indices[jj] == jcol ) {
         *pValue = entries[jj] ;
         return(jj) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jj, nent, ncol, offset, *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   if ( jcol < 0 || jcol >= ncol ) {
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
      offset += sizes[ii] ;
   }
   for ( ii = 0, jj = offset ; ii < sizes[jcol] ; ii++, jj++ ) {
      if ( indices[jj] == irow ) {
         *pValue = entries[jj] ;
         return(jj) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent, *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( irow == rowids[ii] && jcol == colids[ii] ) {
         *pValue = entries[ii] ;
         return(ii) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      ii, joff, nent, nrow, offset, *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( irow < 0 || irow >= nrow || sizes[irow] == 0 ) { 
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < irow ; ii++ ) {
      offset += sizes[ii] ;
   }
   if ( 0 <= (joff = jcol - firstlocs[irow]) && joff < sizes[irow] ) {
      offset += joff ;
      *pValue = entries[offset] ;
      return(offset) ;
   }
   return(-1) ;       
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      ii, ioff, nent, ncol, offset, *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   if ( jcol < 0 || jcol >= ncol || sizes[jcol] == 0 ) { 
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
      offset += sizes[ii] ;
   }
   if ( 0 <= (ioff = irow - firstlocs[jcol]) && ioff < sizes[jcol] ) {
      offset += ioff ;
      *pValue = entries[offset] ;
      return(offset) ;
   }
   return(-1) ;       
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      ncol ;

   if ( irow < 0 || jcol < 0 || irow != jcol ) {
      return(-1) ;
   }
   SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   *pValue = entries[irow] ;
   return(irow) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow < 0 || jcol < 0 ) {
      return(-1) ;
   }
   if ( irow > jcol ) {
      ii   = irow ;
      irow = jcol ;
      jcol = ii   ;
   }
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      size = m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
         if ( jrow == irow ) {
            if ( jcol - irow > m - ii - 1 ) {
               return(-1) ;
            } else {
               kk += jcol - irow ;
               *pValue = entries[kk] ;
               return(kk) ;
            }
         } else {
            kk += size-- ;
         }
      }
   }
   return(kk) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   sign ;
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow < 0 || jcol < 0 ) {
      return(-1) ;
   }
   if ( irow > jcol ) {
      ii   = irow ;
      irow = jcol ;
      jcol = ii   ;
      sign = -1.0 ;
   } else {
      sign = 1.0 ;
   }
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      size = m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
         if ( jrow == irow ) {
            if ( jcol - irow > m - ii - 1 ) {
               return(-1) ;
            } else {
               kk += jcol - irow ;
               *pValue = entries[kk] ;
               return(kk) ;
            }
         } else {
            kk += size-- ;
         }
      }
   }
   return(kk) ;
   } break ;
default :
   fprintf(stderr, 
           "\n fatal error in SubMtx_realEntry(%p,%d,%d,%p)"
           "\n bad mode %d", mtx, irow, jcol, pValue, mtx->mode) ;
   exit(-1) ;
   break ;
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to find matrix entry (irow,jcol) if present.

   return value --
     if entry (irow,jcol) is not present then
        *pReal and *pImag are 0.0
        return value is -1
     else entry (irow,jcol) is present then
        (*pReal,*pImag) is the matrix entry
        return value is offset into entries array 
     endif

   created -- 98may01, cca
   -------------------------------------------------------
*/
int
SubMtx_complexEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || irow < 0 || irow >= mtx->nrow || jcol < 0 
   || jcol >= mtx->ncol || pReal == NULL || pImag == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad input\n", mtx, irow, jcol, pReal, pImag) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtx) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, pReal, pImag, mtx->type) ;
   exit(-1) ;
}
*pReal = *pImag = 0 ;
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, ncol, nrow, offset ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( irow < 0 || irow >= nrow || jcol < 0 || jcol >= ncol ) {
      return(-1) ;
   }
   offset = irow*inc1 + jcol*inc2 ;
   *pReal = entries[2*offset] ;
   *pImag = entries[2*offset+1] ;
   return(offset) ;
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jj, nent, nrow, offset, *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   if ( irow < 0 || irow >= nrow ) {
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < irow ; ii++ ) {
      offset += sizes[ii] ;
   }
   for ( ii = 0, jj = offset ; ii < sizes[irow] ; ii++, jj++ ) {
      if ( indices[jj] == jcol ) {
         *pReal = entries[2*jj] ;
         *pImag = entries[2*jj+1] ;
         return(jj) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jj, nent, ncol, offset, *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   if ( jcol < 0 || jcol >= ncol ) {
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
      offset += sizes[ii] ;
   }
   for ( ii = 0, jj = offset ; ii < sizes[jcol] ; ii++, jj++ ) {
      if ( indices[jj] == irow ) {
         *pReal = entries[2*jj] ;
         *pImag = entries[2*jj+1] ;
         return(jj) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent, *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( irow == rowids[ii] && jcol == colids[ii] ) {
         *pReal = entries[2*ii] ;
         *pImag = entries[2*ii+1] ;
         return(ii) ;
      }
   }
   return(-1) ;
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      ii, joff, nent, nrow, offset, *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( irow < 0 || irow >= nrow || sizes[irow] == 0 ) { 
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < irow ; ii++ ) {
      offset += sizes[ii] ;
   }
   if ( 0 <= (joff = jcol - firstlocs[irow]) && joff < sizes[irow] ) {
      offset += joff ;
      *pReal = entries[2*offset] ;
      *pImag = entries[2*offset+1] ;
      return(offset) ;
   }
   return(-1) ;       
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      ii, ioff, nent, ncol, offset, *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   if ( jcol < 0 || jcol >= ncol || sizes[jcol] == 0 ) { 
      return(-1) ;
   }
   for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
      offset += sizes[ii] ;
   }
   if ( 0 <= (ioff = irow - firstlocs[jcol]) && ioff < sizes[jcol] ) {
      offset += ioff ;
      *pReal = entries[2*offset] ;
      *pImag = entries[2*offset+1] ;
      return(offset) ;
   }
   return(-1) ;       
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      ncol ;

   if ( irow < 0 || jcol < 0 || irow != jcol ) {
      return(-1) ;
   }
   SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   *pReal = entries[2*irow] ;
   *pImag = entries[2*irow+1] ;
   return(irow) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow < 0 || jcol < 0 ) {
      return(-1) ;
   }
   if ( irow > jcol ) {
      ii   = irow ;
      irow = jcol ;
      jcol = ii   ;
   }
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      size = m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
         if ( jrow == irow ) {
            if ( jcol - irow > m - ii - 1 ) {
               return(-1) ;
            } else {
               kk += jcol - irow ;
               *pReal = entries[2*kk] ;
               *pImag = entries[2*kk+1] ;
               return(kk) ;
            }
         } else {
            kk += size-- ;
         }
      }
   }
   return(kk) ;
   } break ;
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   sign ;
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow < 0 || jcol < 0 ) {
      return(-1) ;
   }
   if ( irow > jcol ) {
      ii   = irow ;
      irow = jcol ;
      jcol = ii   ;
      sign = -1.0 ;
   } else {
      sign = 1.0 ;
   }
   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   if ( irow >= ncol || jcol >= ncol ) { 
      return(-1) ;
   }
   for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
      size = m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
         if ( jrow == irow ) {
            if ( jcol - irow > m - ii - 1 ) {
               return(-1) ;
            } else {
               kk += jcol - irow ;
               *pReal = entries[2*kk] ;
               *pImag = sign*entries[2*kk+1] ;
               return(kk) ;
            }
         } else {
            kk += size-- ;
         }
      }
   }
   return(kk) ;
   } break ;
default :
   fprintf(stderr, 
           "\n fatal error in SubMtx_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad mode %d", mtx, irow, jcol, pReal, pImag, mtx->mode) ;
   exit(-1) ;
   break ;
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to return a pointer to the location of
               matrix entry (irow,jcol) if present.

   if entry (irow,jcol) is not present then
      *ppValue is NULL
   else entry (irow,jcol) is present then
      *ppValue is the location of the matrix entry
   endif

   created -- 98may01, cca
   -------------------------------------------------
*/
void
SubMtx_locationOfRealEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   **ppValue
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || irow < 0 || irow >= mtx->nrow || jcol < 0 
   || jcol >= mtx->ncol || ppValue == NULL ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfRealEntry(%p,%d,%d,%p)"
           "\n bad input\n", mtx, irow, jcol, ppValue) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_REAL(mtx) ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfRealEntry(%p,%d,%d,%p)"
           "\n bad type %d, must be SPOOLES_REAL\n", 
           mtx, irow, jcol, ppValue, mtx->type) ;
   exit(-1) ;
}
*ppValue = NULL ;
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, ncol, nrow, offset ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( irow >= 0 && irow < nrow && jcol >= 0 && jcol < ncol ) {
      offset = irow*inc1 + jcol*inc2 ;
      *ppValue = entries + offset ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jj, nent, nrow, offset, *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries);
   if ( irow >= 0 && irow < nrow ) {
      for ( ii = offset = 0 ; ii < irow ; ii++ ) {
         offset += sizes[ii] ;
      }
      for ( ii = 0, jj = offset ; ii < sizes[irow] ; ii++, jj++ ) {
         if ( indices[jj] == jcol ) {
            *ppValue = entries + jj ;
            break ;
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jj, nent, ncol, offset, *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   if ( jcol >= 0 && jcol < ncol ) {
      for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
         offset += sizes[ii] ;
      }
      for ( ii = 0, jj = offset ; ii < sizes[jcol] ; ii++, jj++ ) {
         if ( indices[jj] == irow ) {
            *ppValue = entries + jj ;
            break ;
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent, *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( irow == rowids[ii] && jcol == colids[ii] ) {
         *ppValue = entries + ii ;
         break ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      ii, joff, nent, nrow, offset, *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( irow >= 0 && irow < nrow && sizes[irow] != 0 ) { 
      for ( ii = offset = 0 ; ii < irow ; ii++ ) {
         offset += sizes[ii] ;
      }
      if ( 0 <= (joff = jcol - firstlocs[irow]) 
           && joff < sizes[irow] ) {
         offset += joff ;
         *ppValue = entries + offset ;
         break ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      ii, ioff, nent, ncol, offset, *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   if ( jcol >= 0 && jcol < ncol && sizes[jcol] != 0 ) { 
      for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
         offset += sizes[jcol] ;
      }
      if (  0 <= (ioff = irow - firstlocs[jcol]) 
         && ioff < sizes[jcol] ) {
         offset += ioff ;
         *ppValue = entries + offset ;
         break ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      ncol ;

   if ( irow >= 0 && jcol >= 0 && irow == jcol ) {
      SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
      if ( irow < ncol && jcol < ncol ) { 
         *ppValue = entries + irow ;
      }
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow >= 0 && jcol >= 0 ) {
      SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, 
                               &pivotsizes, &entries) ;
      if ( irow < ncol && jcol < ncol ) { 
         for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
            size = m = pivotsizes[ipivot] ;
            for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
               if ( jrow == irow ) {
                  if ( jrow - irow > m - ii ) {
                     kk = -1 ;
                  } else {
                     kk += jrow - irow ;
                  }
               } else {
                  kk += size-- ;
               }
            }
         }
         if ( kk != -1 ) {
            *ppValue = entries + kk ;
         }
      }
   }
   } break ;
default :
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfRealEntry(%p,%d,%d,%p)"
       "\n bad mode %d", mtx, irow, jcol, ppValue, mtx->mode) ;
   exit(-1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to return a pointer to the location of
               matrix entry (irow,jcol) if present.

   if entry (irow,jcol) is not present then
      (*ppReal,*ppImag) is (NULL,NULL)
   else entry (irow,jcol) is present then
      (*ppReal,*ppImag) is the location of the matrix entry
   endif

   created -- 98may01, cca
   --------------------------------------------------------
*/
void
SubMtx_locationOfComplexEntry (
   SubMtx   *mtx,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL || irow < 0 || irow >= mtx->nrow || jcol < 0 
   || jcol >= mtx->ncol || ppReal == NULL || ppImag == NULL ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad input\n", mtx, irow, jcol, ppReal, ppImag) ;
   exit(-1) ;
}
if ( ! SUBMTX_IS_COMPLEX(mtx) ) {
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfComplexEntry(%p,%d,%d,%p,%p)"
           "\n bad type %d, must be SPOOLES_COMPLEX\n", 
           mtx, irow, jcol, ppReal, ppImag, mtx->type) ;
   exit(-1) ;
}
*ppReal = NULL ;
*ppImag = NULL ;
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS : {
   double   *entries ;
   int      inc1, inc2, ncol, nrow, offset ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   if ( irow >= 0 && irow < nrow && jcol >= 0 && jcol < ncol ) {
      offset = irow*inc1 + jcol*inc2 ;
      *ppReal = entries + 2*offset ;
      *ppImag = entries + 2*offset  + 1 ;
   }
   } break ;
case SUBMTX_SPARSE_ROWS : {
   double   *entries ;
   int      ii, jj, nent, nrow, offset, *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries);
   if ( irow >= 0 && irow < nrow ) {
      for ( ii = offset = 0 ; ii < irow ; ii++ ) {
         offset += sizes[ii] ;
      }
      for ( ii = 0, jj = offset ; ii < sizes[irow] ; ii++, jj++ ) {
         if ( indices[jj] == jcol ) {
            *ppReal = entries + 2*jj ;
            *ppImag = entries + 2*jj  + 1 ;
            break ;
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_COLUMNS : {
   double   *entries ;
   int      ii, jj, nent, ncol, offset, *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   if ( jcol >= 0 && jcol < ncol ) {
      for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
         offset += sizes[ii] ;
      }
      for ( ii = 0, jj = offset ; ii < sizes[jcol] ; ii++, jj++ ) {
         if ( indices[jj] == irow ) {
            *ppReal = entries + 2*jj ;
            *ppImag = entries + 2*jj  + 1 ;
            break ;
         }
      }
   }
   } break ;
case SUBMTX_SPARSE_TRIPLES : {
   double   *entries ;
   int      ii, nent, *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      if ( irow == rowids[ii] && jcol == colids[ii] ) {
         *ppReal = entries + 2*ii ;
         *ppImag = entries + 2*ii  + 1 ;
         break ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBROWS : {
   double   *entries ;
   int      ii, joff, nent, nrow, offset, *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( irow >= 0 && irow < nrow && sizes[irow] != 0 ) { 
      for ( ii = offset = 0 ; ii < irow ; ii++ ) {
         offset += sizes[ii] ;
      }
      if ( 0 <= (joff = jcol - firstlocs[irow]) 
           && joff < sizes[irow] ) {
         offset += joff ;
         *ppReal = entries + 2*offset ;
         *ppImag = entries + 2*offset  + 1 ;
         break ;
      }
   }
   } break ;
case SUBMTX_DENSE_SUBCOLUMNS : {
   double   *entries ;
   int      ii, ioff, nent, ncol, offset, *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   if ( jcol >= 0 && jcol < ncol && sizes[jcol] != 0 ) { 
      for ( ii = offset = 0 ; ii < jcol ; ii++ ) {
         offset += sizes[jcol] ;
      }
      if (  0 <= (ioff = irow - firstlocs[jcol]) 
         && ioff < sizes[jcol] ) {
         offset += ioff ;
         *ppReal = entries + 2*offset ;
         *ppImag = entries + 2*offset  + 1 ;
         break ;
      }
   }
   } break ;
case SUBMTX_DIAGONAL : {
   double   *entries ;
   int      ncol ;

   if ( irow >= 0 && jcol >= 0 && irow == jcol ) {
      SubMtx_diagonalInfo(mtx, &ncol, &entries) ;
      if ( irow < ncol && jcol < ncol ) { 
         *ppReal = entries + 2*irow ;
         *ppImag = entries + 2*irow  + 1 ;
      }
   }
   } break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM : {
   double   *entries ;
   int      ii, ipivot, jrow, kk, m, ncol, nent, size ;
   int      *pivotsizes ;

   if ( irow >= 0 && jcol >= 0 ) {
      SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, 
                               &pivotsizes, &entries) ;
      if ( irow < ncol && jcol < ncol ) { 
         for ( jrow = ipivot = kk = 0 ; jrow <= irow ; ipivot++ ) {
            size = m = pivotsizes[ipivot] ;
            for ( ii = 0 ; ii < m ; ii++, jrow++ ) {
               if ( jrow == irow ) {
                  if ( jrow - irow > m - ii ) {
                     kk = -1 ;
                  } else {
                     kk += jrow - irow ;
                  }
               } else {
                  kk += size-- ;
               }
            }
         }
         if ( kk != -1 ) {
            *ppReal = entries + 2*kk ;
            *ppImag = entries + 2*kk  + 1 ;
         }
      }
   }
   } break ;
default :
   fprintf(stderr, 
       "\n fatal error in SubMtx_locationOfComplexEntry(%p,%d,%d,%p,%p)"
       "\n bad mode %d", mtx, irow, jcol, ppReal, ppImag, mtx->mode) ;
   exit(-1) ;
   break ;
}
return ; }

/*--------------------------------------------------------------------*/
