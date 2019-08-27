/*  instance.c  */

#include "../DenseMtx.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return the row id of the object

   created -- 98may02, cca
   -------------------------------
*/
int
DenseMtx_rowid (
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_rowid(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->rowid) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   return the column id of the object

   created -- 98may02, cca
   ----------------------------------
*/
int
DenseMtx_colid (
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_colid(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->colid) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   fill *pnrow with nrow and *pncol with ncol

   created -- 98may02, cca
   ------------------------------------------
*/
void
DenseMtx_dimensions (
   DenseMtx   *mtx,
   int         *pnrow,
   int         *pncol
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pnrow == NULL || pncol == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_dimensions(%p,%p,%p)"
           "\n bad input\n", mtx, pnrow, pncol) ;
   exit(-1) ;
}
*pnrow = mtx->nrow ;
*pncol = mtx->ncol ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------
   return the row increment

   created -- 98may02, cca
   ------------------------
*/
int
DenseMtx_rowIncrement (
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_rowIncrement(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->inc1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   return the column increment

   created -- 98may02, cca
   ---------------------------
*/
int
DenseMtx_columnIncrement (
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_columnIncrement(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->inc2) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   fill *pnrow with nrow, *prowind with rowind

   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_rowIndices (
   DenseMtx   *mtx,
   int         *pnrow,
   int         **prowind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pnrow == NULL || prowind == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_rowIndices(%p,%p,%p)"
           "\n bad input\n", mtx, pnrow, prowind) ;
   exit(-1) ;
}
*pnrow   = mtx->nrow   ;
*prowind = mtx->rowind ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   fill *pncol with ncol, *pcolind with colind

   created -- 98may02, cca
   -------------------------------------------
*/
void
DenseMtx_columnIndices (
   DenseMtx   *mtx,
   int         *pncol,
   int         **pcolind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pncol == NULL || pcolind == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_columnIndices(%p,%p,%p)"
           "\n bad input\n", mtx, pncol, pcolind) ;
   exit(-1) ;
}
*pncol   = mtx->ncol   ;
*pcolind = mtx->colind ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------
   return a pointer to the entries

   created -- 98may02, cca
   -------------------------------
*/
double *
DenseMtx_entries(
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_entries(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(mtx->entries) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   return a pointer to the workspace

   created -- 98may02, cca
   ---------------------------------
*/
void *
DenseMtx_workspace(
   DenseMtx   *mtx
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_workspace(%p)"
           "\n bad input\n", mtx) ;
   exit(-1) ;
}
return(DV_entries(&mtx->wrkDV)) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   fill *pValue with the entry in (irow,jcol)

   created -- 98jun05, cca
   ------------------------------------------
*/
void
DenseMtx_realEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     *pValue
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pValue == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_realEntry()"
           "\n mtx or pValue is NULL\n") ;
   exit(-1) ;
}
if ( mtx->type != SPOOLES_REAL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_realEntry()"
           "\n mtx type must be SPOOLES_REAL\n") ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_realEntry()"
           "\n irow = %d, mtx->nrow = %d input\n", irow, mtx->nrow) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_realEntry()"
           "\n jcol = %d, mtx->ncol = %d input\n", jcol, mtx->ncol) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_realEntry()"
           "\n mtx->entries is NULL \n") ;
   exit(-1) ;
}
*pValue = mtx->entries[irow*mtx->inc1 + jcol*mtx->inc2] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   fill *pReal and *pImag with the entry in (irow,jcol)

   created -- 98jun05, cca
   ----------------------------------------------------
*/
void
DenseMtx_complexEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     *pReal,
   double     *pImag
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || pReal == NULL || pImag == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_complexEntry()"
           "\n mtxm pReal or pImag is NULL\n") ;
   exit(-1) ;
}
if ( mtx->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in DenseMtx_complexEntry()"
           "\n mtx type must be SPOOLES_COMPLEX\n") ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_complexEntry()"
           "\n irow = %d, mtx->nrow = %d input\n", irow, mtx->nrow) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_complexEntry()"
           "\n jcol = %d, mtx->ncol = %d input\n", jcol, mtx->ncol) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_complexEntry()"
           "\n mtx->entries is NULL \n") ;
   exit(-1) ;
}
loc = 2*(irow*mtx->inc1 + jcol*mtx->inc2) ;
*pReal = mtx->entries[loc]   ;
*pImag = mtx->entries[loc+1] ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   set entry (irow,jcol) to value

   created -- 98jun05, cca
   ------------------------------
*/
void
DenseMtx_setRealEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     value
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setRealEntry()"
           "\n mtx is NULL\n") ;
   exit(-1) ;
}
if ( mtx->type != SPOOLES_REAL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setRealEntry()"
           "\n mtx type must be SPOOLES_REAL\n") ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setRealEntry()"
           "\n irow = %d, mtx->nrow = %d input\n", irow, mtx->nrow) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setRealEntry()"
           "\n jcol = %d, mtx->ncol = %d input\n", jcol, mtx->ncol) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setRealEntry()"
           "\n mtx->entries is NULL \n") ;
   exit(-1) ;
}
mtx->entries[irow*mtx->inc1 + jcol*mtx->inc2] = value ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   set entry (irow,jcol) to (real,imag)

   created -- 98jun05, cca
   ------------------------------------
*/
void
DenseMtx_setComplexEntry (
   DenseMtx   *mtx,
   int        irow,
   int        jcol,
   double     real,
   double     imag
) {
int   loc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setComplexEntry()"
           "\n mtx is NULL\n") ;
   exit(-1) ;
}
if ( mtx->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setComplexEntry()"
           "\n mtx type must be SPOOLES_COMPLEX\n") ;
   exit(-1) ;
}
if ( irow < 0 || irow >= mtx->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setComplexEntry()"
           "\n irow = %d, mtx->nrow = %d input\n", irow, mtx->nrow) ;
   exit(-1) ;
}
if ( jcol < 0 || jcol >= mtx->ncol ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setComplexEntry()"
           "\n jcol = %d, mtx->ncol = %d input\n", jcol, mtx->ncol) ;
   exit(-1) ;
}
if ( mtx->entries == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_setComplexEntry()"
           "\n mtx->entries is NULL \n") ;
   exit(-1) ;
}
loc = 2*(irow*mtx->inc1 + jcol*mtx->inc2) ;
mtx->entries[loc]   = real ;
mtx->entries[loc+1] = imag ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- fill *prowent with the base address of row irow

   return values --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- invalid type for mtx
     -3 -- irow is invalid
     -4 -- prowent is NULL

   created -- 98nov11, cca
   ----------------------------------------------------------
*/
int
DenseMtx_row (
   DenseMtx   *mtx,
   int        irow,
   double     **prowent
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_row()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
if ( mtx->type != SPOOLES_REAL && mtx->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n error in DenseMtx_row()"
           "\n invalid type %d\n", mtx->type) ;
   return(-2) ;
}
if ( irow < 0 || irow >= mtx->nrow ) {
   fprintf(stderr, "\n error in DenseMtx_row()"
           "\n %d rows, irow = %d\n", mtx->nrow, irow) ;
   return(-3) ;
}
if ( prowent == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_row()"
           "\n prowent is NULL\n") ;
   return(-4) ;
}
if ( mtx->type == SPOOLES_REAL ) {
   *prowent = mtx->entries + irow*mtx->inc1 ;
} else {
   *prowent = mtx->entries + 2*irow*mtx->inc1 ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- fill *pcolent with the base address of column jcol

   return values --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- invalid type for mtx
     -3 -- jcol is invalid
     -4 -- pcolent is NULL

   created -- 98nov11, cca
   -------------------------------------------------------------
*/
int
DenseMtx_column (
   DenseMtx   *mtx,
   int        jcol,
   double     **pcolent
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_column()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
if ( mtx->type != SPOOLES_REAL && mtx->type != SPOOLES_COMPLEX ) {
   fprintf(stderr, "\n error in DenseMtx_column()"
           "\n invalid type %d\n", mtx->type) ;
   return(-2) ;
}
if ( jcol < 0 || jcol >= mtx->ncol ) {
   fprintf(stderr, "\n error in DenseMtx_column()"
           "\n %d columns, jcol = %d\n", mtx->ncol, jcol) ;
   return(-3) ;
}
if ( pcolent == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_column()"
           "\n pcolent is NULL\n") ;
   return(-4) ;
}
if ( mtx->type == SPOOLES_REAL ) {
   *pcolent = mtx->entries + jcol*mtx->inc2 ;
} else {
   *pcolent = mtx->entries + 2*jcol*mtx->inc2 ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
