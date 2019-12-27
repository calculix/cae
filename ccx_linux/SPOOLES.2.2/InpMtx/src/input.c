/*  input.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   prepare to add more entries, 
   sort/compress and possible resize

   created -- 98apr25, cca
   ---------------------------------
*/
static void
prepareToAddNewEntries (
   InpMtx   *inpmtx,
   int      nnewent
) {
if ( inpmtx->nent + nnewent > inpmtx->maxnent ) {
/*
   -----------------------------------
   vectors are full, sort and compress
   -----------------------------------
*/
   InpMtx_sortAndCompress(inpmtx) ;
   inpmtx->storageMode = INPMTX_SORTED ;
}
if ( inpmtx->nent + nnewent > inpmtx->maxnent ) {
/*
   ------------------------------
   vectors are still full, resize
   ------------------------------
*/
   int   newmaxnent = inpmtx->maxnent * inpmtx->resizeMultiple ;
   if ( newmaxnent < inpmtx->nent + nnewent ) {
      newmaxnent = inpmtx->nent + nnewent ;
   }
   InpMtx_setMaxnent(inpmtx, newmaxnent) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   input a single entry in the matrix

   created -- 98jan28, cca
   ----------------------------------
*/
static void
inputEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col,
   double    real,
   double    imag
) {
int   nent ;
int   *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, 1) ;
nent  = inpmtx->nent ;
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   ivec1[nent] = row ;
   ivec2[nent] = col ;
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   ivec1[nent] = col ;
   ivec2[nent] = row ;
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   if ( row <= col ) {
      ivec1[nent] = row ;
      ivec2[nent] = col - row ;
   } else {
      ivec1[nent] = col ;
      ivec2[nent] = col - row ;
   }
}
IV_setSize(&inpmtx->ivec1IV, nent + 1) ;
IV_setSize(&inpmtx->ivec2IV, nent + 1) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   dvec[nent] = real ;
   DV_setSize(&inpmtx->dvecDV,  nent + 1) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   dvec[2*nent]   = real  ;
   dvec[2*nent+1] = imag  ;
   DV_setSize(&inpmtx->dvecDV,  2*(nent + 1)) ;
}
inpmtx->nent++ ;
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   input a single entry in the matrix

   created -- 98jan28, cca
   ----------------------------------
*/
void
InpMtx_inputEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col
) {
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL || row < 0 || col < 0 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputEntry(%p,%d,%d)"
           "\n bad input\n", inpmtx, row, col) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputEntry(%p,%d,%d)"
           "\n bad coordType = %d\n", inpmtx, row, col, 
           inpmtx->coordType) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputEntry(%p,%d,%d)"
           "\n inputMode is not INPMTX_INDICES_ONLY\n",
           inpmtx, row, col) ;
   exit(-1) ;
}
inputEntry(inpmtx, row, col, 0.0, 0.0) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   input a single real entry in the matrix

   created -- 98jan28, cca
   ---------------------------------------
*/
void
InpMtx_inputRealEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col,
   double    value
) {
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL || row < 0 || col < 0 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealEntry(%p,%d,%d,%e)"
           "\n bad inputReal\n", inpmtx, row, col, value) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealEntry(%p,%d,%d,%e)"
           "\n bad coordType = %d\n", inpmtx, row, col, value, 
           inpmtx->coordType) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealEntry(%p,%d,%d,%e)"
           "\n inputMode is not SPOOLES_REAL\n",
           inpmtx, row, col, value) ;
   exit(-1) ;
}
inputEntry(inpmtx, row, col, value, 0.0) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   input a single complex entry in the matrix

   created -- 98jan28, cca
   ------------------------------------------
*/
void
InpMtx_inputComplexEntry (
   InpMtx   *inpmtx,
   int       row,
   int       col,
   double    real,
   double    imag
) {
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL || row < 0 || col < 0 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexEntry(%p,%d,%d,%e,%e)"
           "\n bad inputComplex\n", inpmtx, row, col, real, imag) ;
   exit(-1) ;
}
if ( !(   INPMTX_IS_BY_ROWS(inpmtx)
       || INPMTX_IS_BY_COLUMNS(inpmtx)
       || INPMTX_IS_BY_CHEVRONS(inpmtx) ) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexEntry(%p,%d,%d,%e,%e)"
           "\n bad coordType = %d\n", inpmtx, row, col, real, imag,
           inpmtx->coordType) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexEntry(%p,%d,%d,%e,%e)"
           "\n inputMode is not SPOOLES_COMPLEX\n",
           inpmtx, row, col, real, imag) ;
   exit(-1) ;
}
inputEntry(inpmtx, row, col, real, imag) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   input a row in the matrix

   created -- 98jan28, cca
   ---------------------------------
*/
static void
inputRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[],
   double    rowent[]
) {
int      col, ii, jj, nent ;
int      *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, rowsize) ;
nent  = inpmtx->nent ; 
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) { /* row coordinates */
   IVfill(rowsize, ivec1 + nent, row) ;
   IVcopy(rowsize, ivec2 + nent, rowind) ;
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) { /* column coordinates */
   IVfill(rowsize, ivec2 + nent, row) ;
   IVcopy(rowsize, ivec1 + nent, rowind) ;
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) { /* chevron coordinates */
   for ( ii = 0, jj = nent ; ii < rowsize ; ii++, jj++ ) {
      col = rowind[ii] ;
      ivec1[ii] = (row <= col) ? row : col ;
      ivec2[ii] = col - row ;
   }
}
IV_setSize(&inpmtx->ivec1IV, nent + rowsize) ;
IV_setSize(&inpmtx->ivec2IV, nent + rowsize) ;
/*
   -----------------
   input the entries
   -----------------
*/
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double  *dvec = DV_entries(&inpmtx->dvecDV) ;
   DVcopy(rowsize, dvec + nent, rowent) ;
   DV_setSize(&inpmtx->dvecDV, nent + rowsize) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double  *dvec = DV_entries(&inpmtx->dvecDV) ;
   ZVcopy(rowsize, dvec + 2*nent, rowent) ;
   DV_setSize(&inpmtx->dvecDV,  2*(nent + rowsize)) ;
}
inpmtx->storageMode = INPMTX_RAW_DATA ;
inpmtx->nent += rowsize ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   input a real row in the matrix

   created -- 98jan28, cca
   ------------------------------
*/
void
InpMtx_inputRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || row < 0 || rowsize < 0 || rowind == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRow(%p,%d,%d,%p)"
           "\n bad input\n", inpmtx, row, rowsize, rowind) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRow(%p,%d,%d,%p)"
           "\n inputMode is not INPMTX_INDICES_ONLY\n",
           inpmtx, row, rowsize, rowind) ;
   exit(-1) ;
}
inputRow(inpmtx, row, rowsize, rowind, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------
   input a real row in the matrix

   created -- 98jan28, cca
   ------------------------------
*/
void
InpMtx_inputRealRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[],
   double    rowent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || row < 0 || rowsize < 0 
   || rowind == NULL || rowent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealRow(%p,%d,%d,%p,%p)"
           "\n bad input\n", inpmtx, row, rowsize, rowind, rowent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealRow(%p,%d,%d,%p,%p)"
           "\n inputMode is not SPOOLES_REAL\n",
           inpmtx, row, rowsize, rowind, rowent) ;
   exit(-1) ;
}
inputRow(inpmtx, row, rowsize, rowind, rowent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   input a complex row in the matrix

   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_inputComplexRow (
   InpMtx   *inpmtx,
   int       row,
   int       rowsize,
   int       rowind[],
   double    rowent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || row < 0 || rowsize < 0 
   || rowind == NULL || rowent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexRow(%p,%d,%d,%p,%p)"
           "\n bad input\n", inpmtx, row, rowsize, rowind, rowent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexRow(%p,%d,%d,%p,%p)"
           "\n inputMode is not SPOOLES_COMPLEX\n",
           inpmtx, row, rowsize, rowind, rowent) ;
   exit(-1) ;
}
inputRow(inpmtx, row, rowsize, rowind, rowent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   input a complex column in the matrix

   created -- 98jan28, cca
   ------------------------------------
*/
static void
inputColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[],
   double    colent[]
) {
int      ii, jj, nent, row ;
int      *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, colsize) ;
nent  = inpmtx->nent ;
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   IVcopy(colsize, ivec1 + nent, colind) ;
   IVfill(colsize, ivec2 + nent, col) ;
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   IVfill(colsize, ivec1 + nent, col) ;
   IVcopy(colsize, ivec2 + nent, colind) ;
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( ii = 0, jj = nent ; ii < colsize ; ii++, jj++ ) {
      row = colind[jj] ;
      ivec1[jj] = (row <= col) ? row : col ;
      ivec2[jj] = col - row ;
   }
}
IV_setSize(&inpmtx->ivec1IV, nent + colsize) ;
IV_setSize(&inpmtx->ivec2IV, nent + colsize) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double *dvec = DV_entries(&inpmtx->dvecDV) + nent ;
   DVcopy(colsize, dvec, colent) ;
   DV_setSize(&inpmtx->dvecDV,  nent + colsize) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double *dvec = DV_entries(&inpmtx->dvecDV) + 2*nent ;
   ZVcopy(colsize, dvec, colent) ;
   DV_setSize(&inpmtx->dvecDV,  2*(nent + colsize)) ;
}
inpmtx->nent = nent + colsize ;
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   input a column in the matrix

   created -- 98jan28, cca
   ----------------------------
*/
void
InpMtx_inputColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || col < 0 || colsize < 0 || colind == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealColumn(%p,%d,%d,%p)"
           "\n bad input\n", inpmtx, col, colsize, colind) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputColumn(%p,%d,%d,%p)"
           "\n inputMode must be INPMTX_INDICES_ONLY\n",
           inpmtx, col, colsize, colind) ;
   exit(-1) ;
}
inputColumn(inpmtx, col, colsize, colind, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   input a real column in the matrix

   created -- 98jan28, cca
   ---------------------------------
*/
void
InpMtx_inputRealColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[],
   double    colent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || col < 0 || colsize < 0 
   || colind == NULL || colent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealColumn(%p,%d,%d,%p,%p)"
           "\n bad input\n", inpmtx, col, colsize, colind, colent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealColumn(%p,%d,%d,%p,%p)"
           "\n inputMode must be SPOOLES_REAL\n",
           inpmtx, col, colsize, colind, colent) ;
   exit(-1) ;
}
inputColumn(inpmtx, col, colsize, colind, colent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   input a complex column in the matrix

   created -- 98jan28, cca
   ------------------------------------
*/
void
InpMtx_inputComplexColumn (
   InpMtx   *inpmtx,
   int       col,
   int       colsize,
   int       colind[],
   double    colent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || col < 0 || colsize < 0 
   || colind == NULL || colent == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexColumn(%p,%d,%d,%p,%p)"
           "\n bad input\n", inpmtx, col, colsize, colind, colent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexColumn(%p,%d,%d,%p,%p)"
           "\n inputMode must be SPOOLES_COMPLEX\n",
           inpmtx, col, colsize, colind, colent) ;
   exit(-1) ;
}
inputColumn(inpmtx, col, colsize, colind, colent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   input a chevron in the matrix

   created -- 98jan28, cca
   -----------------------------
*/
static void
inputChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[],
   double    chvent[]
) {
int      col, ii, jj, nent, offset, row ;
int      *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, chvsize) ;
nent  = inpmtx->nent ;
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( ii = 0, jj = nent ; ii < chvsize ; ii++, jj++ ) {
      if ( (offset = chvind[ii]) >= 0 ) {
         row = chv ;
         col = chv + offset ;
      } else {
         col = chv ;
         row = chv - offset ;
      }
      ivec1[jj] = row ;
      ivec2[jj] = col ;
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( ii = 0, jj = nent ; ii < chvsize ; ii++, jj++ ) {
      if ( (offset = chvind[ii]) >= 0 ) {
         row = chv ;
         col = chv + offset ;
      } else {
         col = chv ;
         row = chv - offset ;
      }
      ivec1[jj] = col ;
      ivec2[jj] = row ;
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   IVfill(chvsize, ivec1 + nent, chv) ;
   IVcopy(chvsize, ivec2 + nent, chvind) ;
}
IV_setSize(&inpmtx->ivec1IV, nent + chvsize) ;
IV_setSize(&inpmtx->ivec2IV, nent + chvsize) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) + nent ;
   DVcopy(chvsize, dvec, chvent) ;
   DV_setSize(&inpmtx->dvecDV,  nent + chvsize) ;
} else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) + 2*nent ;
   ZVcopy(chvsize, dvec, chvent) ;
   DV_setSize(&inpmtx->dvecDV,  2*(nent + chvsize)) ;
}
inpmtx->nent += chvsize ;
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   input a chevron in the matrix

   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || chv < 0 || chvsize < 0 || chvind == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputChevron(%p,%d,%d,%p)"
          "\n bad input\n", inpmtx, chv, chvsize, chvind) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputChevron(%p,%d,%d,%p)"
          "\n inputMode must be INPMTX_INDICES_ONLY\n", 
           inpmtx, chv, chvsize, chvind) ;
   exit(-1) ;
}
inputChevron(inpmtx, chv, chvsize, chvind, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   input a chevron in the matrix

   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputRealChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[],
   double    chvent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || chv < 0 || chvsize < 0 
   || chvind == NULL || chvent == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputRealChevron(%p,%d,%d,%p,%p)"
          "\n bad input\n", inpmtx, chv, chvsize, chvind, chvent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputRealChevron(%p,%d,%d,%p,%p)"
          "\n inputMode must be SPOOLES_REAL\n", 
           inpmtx, chv, chvsize, chvind, chvent) ;
   exit(-1) ;
}
inputChevron(inpmtx, chv, chvsize, chvind, chvent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   input a chevron in the matrix

   created -- 98jan28, cca
   -----------------------------
*/
void
InpMtx_inputComplexChevron (
   InpMtx   *inpmtx,
   int       chv,
   int       chvsize,
   int       chvind[],
   double    chvent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || chv < 0 || chvsize < 0 
   || chvind == NULL || chvent == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputComplexChevron(%p,%d,%d,%p,%p)"
          "\n bad input\n", inpmtx, chv, chvsize, chvind, chvent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputComplexChevron(%p,%d,%d,%p,%p)"
          "\n inputMode must be SPOOLES_COMPLEX\n", 
           inpmtx, chv, chvsize, chvind, chvent) ;
   exit(-1) ;
}
inputChevron(inpmtx, chv, chvsize, chvind, chvent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   input a matrix

   created -- 98jan28, cca
   -----------------------
*/
static void
inputMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[],
   double    mtxent[]
) {
int      col, ii, jj, kk, nent, row ;
int      *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, nrow*ncol) ;
nent  = inpmtx->nent ;
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   for ( jj = 0, kk = nent ; jj < ncol ; jj++ ) {
      col = colind[jj] ;
      for ( ii = 0 ; ii < nrow ; ii++, kk++ ) {
         row = rowind[ii] ;
         ivec1[kk] = row ;
         ivec2[kk] = col ;
      }
   }
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   for ( jj = 0, kk = nent ; jj < ncol ; jj++ ) {
      col = colind[jj] ;
      for ( ii = 0 ; ii < nrow ; ii++, kk++ ) {
         row = rowind[ii] ;
         ivec1[kk] = col ;
         ivec2[kk] = row ;
      }
   }
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   for ( jj = 0, kk = nent ; jj < ncol ; jj++ ) {
      col = colind[jj] ;
      for ( ii = 0 ; ii < nrow ; ii++, kk++ ) {
         row = rowind[ii] ;
         if ( row <= col ) {
            ivec1[kk] = row ;
         } else {
            ivec1[kk] = col ;
         }
         ivec2[kk] = col - row ;
      }
   }
}
IV_setSize(&inpmtx->ivec1IV, nent + nrow*ncol) ;
IV_setSize(&inpmtx->ivec2IV, nent + nrow*ncol) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   int      ij ;
   for ( jj = 0, kk = nent ; jj < ncol ; jj++ ) {
      for ( ii = 0 ; ii < nrow ; ii++, kk++ ) {
         ij = ii*rowstride + jj*colstride ;
         dvec[kk] = mtxent[ij] ;
      }
   }
   DV_setSize(&inpmtx->dvecDV, nent + nrow*ncol) ;
} if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   int      ij ;
   for ( jj = 0, kk = nent ; jj < ncol ; jj++ ) {
      for ( ii = 0 ; ii < nrow ; ii++, kk++ ) {
         ij = ii*rowstride + jj*colstride ;
         dvec[2*kk]   = mtxent[2*ij]   ;
         dvec[2*kk+1] = mtxent[2*ij+1] ;
      }
   }
   DV_setSize(&inpmtx->dvecDV,  2*(nent + nrow*ncol)) ;
}
inpmtx->nent += nrow*ncol ;
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   input a matrix

   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || nrow < 0 || ncol < 0 
   || rowstride < 1  || colstride < 1
   || rowind == NULL || colind == NULL ) {
   fprintf(stderr, 
  "\n fatal error in InpMtx_inputMatrix(%p,%d,%d,%d,%d,%p,%p)"
  "\n bad input\n", inpmtx, nrow, ncol, rowstride, colstride, 
        rowind, colind) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
 "\n fatal error in InpMtx_inputMatrix(%p,%d,%d,%d,%d,%p,%p)"
 "\n inputComplexMode must be INPMTX_INDICES_ONLY\n",
        inpmtx, nrow, ncol, rowstride, colstride, rowind, colind) ;
   exit(-1) ;
}
if ( nrow == 0 || ncol == 0 ) {
   return ;
}
inputMatrix(inpmtx, nrow, ncol, rowstride, colstride, 
            rowind, colind, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   input a matrix

   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputRealMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[],
   double    mtxent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || nrow < 0 || ncol < 0 
   || rowstride < 1 || colstride < 1
   || rowind == NULL || colind == NULL || mtxent == NULL ) {
   fprintf(stderr, 
  "\n fatal error in InpMtx_inputRealMatrix(%p,%d,%d,%d,%d,%p,%p,%p)"
  "\n bad input\n", inpmtx, nrow, ncol, rowstride, colstride, 
        rowind, colind, mtxent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
 "\n fatal error in InpMtx_inputRealMatrix(%p,%d,%d,%d,%d,%p,%p,%p)"
 "\n inputMode must be SPOOLES_REAL\n",
        inpmtx, nrow, ncol, rowstride, colstride, 
        rowind, colind, mtxent) ;
   exit(-1) ;
}
if ( nrow == 0 || ncol == 0 ) {
   return ;
}
inputMatrix(inpmtx, nrow, ncol, rowstride, colstride, 
            rowind, colind, mtxent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------
   input a matrix

   created -- 98jan28, cca
   -----------------------
*/
void
InpMtx_inputComplexMatrix (
   InpMtx   *inpmtx,
   int       nrow,
   int       ncol,
   int       rowstride,
   int       colstride,
   int       rowind[],
   int       colind[],
   double    mtxent[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || nrow < 0 || ncol < 0 
   || rowstride < 1 || colstride < 1
   || rowind == NULL || colind == NULL || mtxent == NULL ) {
   fprintf(stderr, 
  "\n fatal error in InpMtx_inputComplexMatrix(%p,%d,%d,%d,%d,%p,%p,%p)"
  "\n bad input\n", inpmtx, nrow, ncol, rowstride, colstride, 
        rowind, colind, mtxent) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
 "\n fatal error in InpMtx_inputComplexMatrix(%p,%d,%d,%d,%d,%p,%p,%p)"
 "\n inputMode must be SPOOLES_COMPLEX\n",
        inpmtx, nrow, ncol, rowstride, colstride, 
        rowind, colind, mtxent) ;
   exit(-1) ;
}
if ( nrow == 0 || ncol == 0 ) {
   return ;
}
inputMatrix(inpmtx, nrow, ncol, rowstride, colstride, 
            rowind, colind, mtxent) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   inputComplex a number of (row,column, entry) triples into the matrix

   created -- 98jan28, cca
   --------------------------------------------------------------------
*/
static void
inputTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[],
   double    entries[]
) {
int      nent ;
int      *ivec1, *ivec2 ;

prepareToAddNewEntries(inpmtx, ntriples) ;
nent  = inpmtx->nent ;
ivec1 = IV_entries(&inpmtx->ivec1IV) ;
ivec2 = IV_entries(&inpmtx->ivec2IV) ;
IVcopy(ntriples, ivec1 + nent, rowids) ;
IVcopy(ntriples, ivec2 + nent, colids) ;
IV_setSize(&inpmtx->ivec1IV, nent + ntriples) ;
IV_setSize(&inpmtx->ivec2IV, nent + ntriples) ;
if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   DVcopy(ntriples, dvec + nent, entries) ;
   DV_setSize(&inpmtx->dvecDV,  nent + ntriples) ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   double   *dvec = DV_entries(&inpmtx->dvecDV) ;
   ZVcopy(ntriples, dvec + 2*nent, entries) ;
   DV_setSize(&inpmtx->dvecDV,  2*(nent + ntriples)) ;
}
inpmtx->nent += ntriples ;
inpmtx->storageMode = INPMTX_RAW_DATA ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[]
) {
/*
   --------------
   check the data
   --------------
*/
if (  inpmtx == NULL || ntriples < 0 
   || rowids == NULL || colids == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputTriples(%p,%d,%p,%p)"
          "\n bad inputComplex\n", 
          inpmtx, ntriples, rowids, colids) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputEntry(%p,%d,%p,%p)"
           "\n coordType must be INPMTX_INDICES_ONLY\n",
           inpmtx, ntriples, rowids, colids) ;
   exit(-1) ;
}
inputTriples(inpmtx, ntriples, rowids, colids, NULL) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputRealTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[],
   double    entries[]
) {
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL || ntriples < 0 || rowids == NULL 
   || colids == NULL || entries == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputRealTriples(%p,%d,%p,%p,%p)"
          "\n bad input\n", 
          inpmtx, ntriples, rowids, colids, entries) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputRealEntry(%p,%d,%p,%p,%p)"
           "\n coordType must be COMPLEX_REAL_ENTRIES\n",
           inpmtx, ntriples, rowids, colids, entries) ;
   exit(-1) ;
}
inputTriples(inpmtx, ntriples, rowids, colids, entries) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   input a number of (row,column, entry) triples into the matrix

   created -- 98jan28, cca
   -------------------------------------------------------------
*/
void
InpMtx_inputComplexTriples (
   InpMtx   *inpmtx,
   int       ntriples,
   int       rowids[],
   int       colids[],
   double    entries[]
) {
/*
   --------------
   check the data
   --------------
*/
if ( inpmtx == NULL || ntriples < 0 || rowids == NULL 
   || colids == NULL || entries == NULL ) {
   fprintf(stderr, 
          "\n fatal error in InpMtx_inputComplexTriples(%p,%d,%p,%p,%p)"
          "\n bad inputComplex\n", 
          inpmtx, ntriples, rowids, colids, entries) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_inputComplexEntry(%p,%d,%p,%p,%p)"
           "\n inputMode must be SPOOLES_COMPLEX\n",
           inpmtx, ntriples, rowids, colids, entries) ;
   exit(-1) ;
}
inputTriples(inpmtx, ntriples, rowids, colids, entries) ;

return ; }

/*--------------------------------------------------------------------*/
