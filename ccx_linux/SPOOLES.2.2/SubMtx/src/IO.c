/*  IO.c  */

#include "../SubMtx.h"

/*--------------------------------------------------------------------*/

static const char *suffixb = ".submtxb" ;
static const char *suffixf = ".submtxf" ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to read in an SubMtx object from a file

   input --

      fn -- filename, must be *.submtxb or *.submtxf

   return value -- 1 if success, 0 if failure

   created -- 98feb15, cca
   --------------------------------------------------
*/
int
SubMtx_readFromFile ( 
   SubMtx   *mtx, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in SubMtx_readFromFile(%p,%s)"
           "\n bad input\n", mtx, fn) ;
   return(0) ;
}
/*
   -------------
   read the file
   -------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "rb")) == NULL ) {
         fprintf(stderr, "\n error in SubMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", mtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SubMtx_readFromBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in SubMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", mtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SubMtx_readFromFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in SubMtx_readFromFile(%p,%s)"
              "\n bad SubMtx file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              mtx, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in SubMtx_readFromFile(%p,%s)"
       "\n bad SubMtx file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       mtx, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to read an SubMtx object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   -------------------------------------------------------
*/
int
SubMtx_readFromFormattedFile ( 
   SubMtx    *mtx, 
   FILE   *fp 
) {
double   *entries ;
int      inc1, inc2, ncol, nent, nrow, rc ;
int      itemp[7] ;
int      *colids, *colind, *firstlocs, *indices, 
         *pivotsizes, *rowids, *rowind, *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   return(0) ;
}
SubMtx_clearData(mtx) ;
/*
   -------------------------------
   read in the seven scalar fields
   -------------------------------
*/
if ( (rc = IVfscanf(fp, 7, itemp)) != 7 ) {
   fprintf(stderr, "\n 1. error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, 7) ;
   return(0) ;
}
/*
   ---------------------
   check the matrix type
   ---------------------
*/
switch ( itemp[0] ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n type = %d not supported\n", mtx, fp, itemp[0]) ;
   return(0) ;
}
/*
   ---------------------
   check the matrix mode
   ---------------------
*/
switch ( itemp[1] ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_DIAGONAL :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   break ;
default :
   fprintf(stderr, "\n error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n mode = %d not supported\n", mtx, fp, itemp[1]) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
SubMtx_init(mtx, itemp[0], itemp[1], itemp[2], 
            itemp[3], itemp[4], itemp[5], itemp[6]) ;
/*
   ---------------------------
   read in the rowind[] vector
   ---------------------------
*/
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( (rc = IVfscanf(fp, nrow, rowind)) != nrow ) {
   fprintf(stderr, "\n 2. error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
   return(0) ;
}
/*
   ---------------------------
   read in the colind[] vector
   ---------------------------
*/
SubMtx_columnIndices(mtx, &ncol, &colind) ;
if ( (rc = IVfscanf(fp, ncol, colind)) != ncol ) {
   fprintf(stderr, "\n 3. error in SubMtx_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
   return(0) ;
}
/*
   ---------------------------------------------------------------
   get pointers and dimensions and read in any integer information
   ---------------------------------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   nent = nrow*ncol ;
   break ;
case SUBMTX_SPARSE_ROWS :
   SubMtx_sparseRowsInfo(mtx, &nrow, &nent, &sizes, &indices, &entries) ;
   if ( (rc = IVfscanf(fp, nrow, sizes)) != nrow ) {
      fprintf(stderr, 
              "\n 5. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   if ( (rc = IVfscanf(fp, nent, indices)) != nent ) {
      fprintf(stderr, 
              "\n 6. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_SPARSE_COLUMNS :
   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   if ( (rc = IVfscanf(fp, ncol, sizes)) != ncol ) {
      fprintf(stderr, 
              "\n 8. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   if ( (rc = IVfscanf(fp, nent, indices)) != nent ) {
      fprintf(stderr, 
              "\n 6. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items9of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_SPARSE_TRIPLES :
   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   if ( (rc = IVfscanf(fp, nent, rowids)) != nent ) {
      fprintf(stderr, 
              "\n 11. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   if ( (rc = IVfscanf(fp, nent, colids)) != nent ) {
      fprintf(stderr, 
              "\n 12. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_DENSE_SUBROWS :
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( (rc = IVfscanf(fp, nrow, firstlocs)) != nrow ) {
      fprintf(stderr, 
              "\n 14. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   if ( (rc = IVfscanf(fp, nrow, sizes)) != nrow ) {
      fprintf(stderr, 
              "\n 15. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( (rc = IVfscanf(fp, ncol, firstlocs)) != ncol ) {
      fprintf(stderr, 
              "\n 14. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   if ( (rc = IVfscanf(fp, ncol, sizes)) != ncol ) {
      fprintf(stderr, 
              "\n 15. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   break ;
case SUBMTX_DIAGONAL :
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   if ( (rc = IVfscanf(fp, nrow, pivotsizes)) != nrow ) {
      fprintf(stderr, 
              "\n 16. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   break ;
}
/*
   --------------------------
   read in the matrix entries
   --------------------------
*/
if ( SUBMTX_IS_REAL(mtx) ) {
   if ( (rc = DVfscanf(fp, nent, entries)) != nent ) {
      fprintf(stderr, 
              "\n 4. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   if ( (rc = DVfscanf(fp, 2*nent, entries)) != 2*nent ) {
      fprintf(stderr, 
              "\n 4. error in SubMtx_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, 2*nent) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read an SubMtx object from a binary file

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   ------------------------------------------------------
*/
int
SubMtx_readFromBinaryFile ( 
   SubMtx   *mtx, 
   FILE     *fp 
) {
double   *entries ;
int      inc1, inc2, ncol, nent, nrow, rc ;
int      itemp[7] ;
int      *colids, *colind, *firstlocs, *indices, 
         *pivotsizes, *rowids, *rowind, *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   return(0) ;
}
SubMtx_clearData(mtx) ;
/*
   -------------------------------
   read in the seven scalar fields
   -------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 7, fp)) != 7 ) {
   fprintf(stderr, "\n 1. error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, 7) ;
   return(0) ;
}
/*
   ---------------------
   check the matrix type
   ---------------------
*/
switch ( itemp[0] ) {
case SPOOLES_REAL :
case SPOOLES_COMPLEX :
   break ;
default :
   fprintf(stderr, "\n error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n type = %d not supported\n", mtx, fp, itemp[0]) ;
   return(0) ;
}
/*
   ---------------------
   check the matrix mode
   ---------------------
*/
switch ( itemp[1] ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
case SUBMTX_SPARSE_ROWS :
case SUBMTX_SPARSE_COLUMNS :
case SUBMTX_DENSE_SUBROWS :
case SUBMTX_DENSE_SUBCOLUMNS :
case SUBMTX_DIAGONAL :
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   break ;
default :
   fprintf(stderr, "\n error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n mode = %d not supported\n", mtx, fp, itemp[1]) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
SubMtx_init(mtx, itemp[0], itemp[1], itemp[2], 
            itemp[3], itemp[4], itemp[5], itemp[6]) ;
/*
   ---------------------------
   read in the rowind[] vector
   ---------------------------
*/
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( (rc = fread((void *) rowind, sizeof(int), nrow, fp)) != nrow) {
   fprintf(stderr, "\n 2. error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
   return(0) ;
}
/*
   ---------------------------
   read in the colind[] vector
   ---------------------------
*/
SubMtx_columnIndices(mtx, &ncol, &colind) ;
if ( (rc = fread((void *) colind, sizeof(int), ncol, fp)) != ncol) {
   fprintf(stderr, "\n 3. error in SubMtx_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
   return(0) ;
}
/*
   ---------------------------------------------------------------
   get pointers and dimensions and read in any integer information
   ---------------------------------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   nent = nrow*ncol ;
   break ;
case SUBMTX_SPARSE_ROWS :
   SubMtx_sparseRowsInfo(mtx, 
                         &nrow, &nent, &sizes, &indices, &entries) ;
   if ( (rc = fread((void *) sizes, sizeof(int), nrow, fp)) != nrow) {
      fprintf(stderr, "\n 5. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   if ( (rc = fread((void *) indices, sizeof(int), nent, fp)) != nent) {
      fprintf(stderr, "\n 6. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_SPARSE_COLUMNS :
   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                            &sizes, &indices, &entries) ;
   if ( (rc = fread((void *) sizes, sizeof(int), ncol, fp)) != ncol) {
      fprintf(stderr, "\n 8. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   if ( (rc = fread((void *) indices, sizeof(int), nent, fp)) != nent) {
      fprintf(stderr, "\n 6. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items9of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_SPARSE_TRIPLES :
   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   if ( (rc = fread((void *) rowids, sizeof(int), nent, fp)) != nent) {
      fprintf(stderr, 
              "\n 11. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   if ( (rc = fread((void *) colids, sizeof(int), nent, fp)) != nent) {
      fprintf(stderr, 
              "\n 12. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
   break ;
case SUBMTX_DENSE_SUBROWS :
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( (rc = fread((void *) firstlocs, sizeof(int), nrow, fp)) 
        != nrow) {
      fprintf(stderr, 
              "\n 14. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   if ( (rc = fread((void *) sizes, sizeof(int), nrow, fp)) 
        != nrow) {
      fprintf(stderr, 
              "\n 15. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                         &firstlocs, &sizes, &entries) ;
   if ( (rc = fread((void *) firstlocs, sizeof(int), ncol, fp)) 
        != ncol) {
      fprintf(stderr, 
              "\n 14. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   if ( (rc = fread((void *) sizes, sizeof(int), ncol, fp)) 
        != ncol) {
      fprintf(stderr, 
              "\n 15. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, ncol) ;
      return(0) ;
   }
   break ;
case SUBMTX_DIAGONAL :
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   if ( (rc = fread((void *) pivotsizes, sizeof(int), nrow, fp)) 
        != nrow) {
      fprintf(stderr, 
              "\n 16. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nrow) ;
      return(0) ;
   }
   break ;
}
/*
   --------------------------
   read in the matrix entries
   --------------------------
*/
if ( SUBMTX_IS_REAL(mtx) ) {
   if ( (rc = fread((void *) entries, sizeof(double), nent, fp)) 
       != nent) {
      fprintf(stderr, "\n 4. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, nent) ;
      return(0) ;
   }
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   if ( (rc = fread((void *) entries, sizeof(double), 2*nent, fp)) 
       != 2*nent) {
      fprintf(stderr, "\n 4. error in SubMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", mtx, fp, rc, 2*nent) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to write an SubMtx object to a file

   input --

      fn -- filename
        *.submtxb -- binary
        *.submtxf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   ----------------------------------------------
*/
int
SubMtx_writeToFile ( 
   SubMtx   *mtx, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeToFile(%p,%s)"
    "\n bad input\n", mtx, fn) ; 
}
/*
   ------------------
   write out the file
   ------------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "wb")) == NULL ) {
         fprintf(stderr, "\n error in SubMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", mtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SubMtx_writeToBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in SubMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", mtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SubMtx_writeToFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in SubMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", mtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SubMtx_writeForHumanEye(mtx, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in SubMtx_writeToFile(%p,%s)"
              "\n unable to open file %s", mtx, fn, fn) ;
      rc = 0 ;
   } else {
      rc = SubMtx_writeForHumanEye(mtx, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to write an SubMtx object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   --------------------------------------------------------
*/
int
SubMtx_writeToFormattedFile ( 
   SubMtx   *mtx, 
   FILE     *fp 
) {
double   *entries ;
int      inc1, inc2, ncol, nent, nrow ;
int      itemp[7] ;
int      *colids, *colind, *firstlocs, *indices, 
         *pivotsizes, *rowids, *rowind, *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeToFormattedFile(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
/*
    ---------------------------
    write out the scalar values
    ---------------------------
*/
itemp[0] = mtx->type  ;
itemp[1] = mtx->mode  ;
itemp[2] = mtx->rowid ;
itemp[3] = mtx->colid ;
itemp[4] = mtx->nrow  ;
itemp[5] = mtx->ncol  ;
itemp[6] = mtx->nent  ;
IVfprintf(fp, 7, itemp) ;
/*
   ------------------------------------
   write out the row and column indices
   ------------------------------------
*/
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
IVfprintf(fp, nrow, rowind) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
IVfprintf(fp, ncol, colind) ;
/*
   ---------------------------------------------------------------------
   get the dimensions and pointers and write out any integer information
   ---------------------------------------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   nent = nrow*ncol ;
   break ;
case SUBMTX_SPARSE_ROWS :
   SubMtx_sparseRowsInfo(mtx, 
                         &nrow, &nent, &sizes, &indices, &entries) ;
   IVfprintf(fp, nrow, sizes) ;
   IVfprintf(fp, nent, indices) ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                            &sizes, &indices, &entries) ;
   IVfprintf(fp, ncol, sizes) ;
   IVfprintf(fp, nent, indices) ;
   break ;
case SUBMTX_SPARSE_TRIPLES :
   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   IVfprintf(fp, nent, rowids) ;
   IVfprintf(fp, nent, colids) ;
   break ;
case SUBMTX_DENSE_SUBROWS :
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                           &firstlocs, &sizes, &entries) ;
   IVfprintf(fp, nrow, firstlocs) ;
   IVfprintf(fp, nrow, sizes) ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                              &firstlocs, &sizes, &entries) ;
   IVfprintf(fp, ncol, firstlocs) ;
   IVfprintf(fp, ncol, sizes) ;
   break ;
case SUBMTX_DIAGONAL :
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   IVfprintf(fp, nrow, pivotsizes) ;
   break ;
}
/*
   ---------------------
   write out the entries
   ---------------------
*/
if ( SUBMTX_IS_REAL(mtx) ) {
   DVfprintf(fp, nent, entries) ;
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   DVfprintf(fp, 2*nent, entries) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an SubMtx object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
SubMtx_writeToBinaryFile ( 
   SubMtx   *mtx, 
   FILE     *fp 
) {
double   *entries ;
int      inc1, inc2, ncol, nent, nrow, rc ;
int      itemp[7] ;
int      *colids, *colind, *firstlocs, *indices, 
         *pivotsizes, *rowids, *rowind, *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeToBinaryFile(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
/*
    ---------------------------
    write out the scalar values
    ---------------------------
*/
itemp[0] = mtx->type  ;
itemp[1] = mtx->mode  ;
itemp[2] = mtx->rowid ;
itemp[3] = mtx->colid ;
itemp[4] = mtx->nrow  ;
itemp[5] = mtx->ncol  ;
itemp[6] = mtx->nent  ;
rc = fwrite((void *) itemp, sizeof(int), 7, fp) ;
/*
   ------------------------------------
   write out the row and column indices
   ------------------------------------
*/
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
rc = fwrite((void *) rowind, sizeof(int), nrow, fp) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
rc = fwrite((void *) colind, sizeof(int), ncol, fp) ;
/*
   ---------------------------------------------------------------------
   get the dimensions and pointers and write out any integer information
   ---------------------------------------------------------------------
*/
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
case SUBMTX_DENSE_COLUMNS :
   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   nent = nrow*ncol ;
   break ;
case SUBMTX_SPARSE_ROWS :
   SubMtx_sparseRowsInfo(mtx, 
                         &nrow, &nent, &sizes, &indices, &entries) ;
   rc = fwrite((void *) sizes,   sizeof(int), nrow, fp) ;
   rc = fwrite((void *) indices, sizeof(int), nent, fp) ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                            &sizes, &indices, &entries) ;
   rc = fwrite((void *) sizes,   sizeof(int), ncol, fp) ;
   rc = fwrite((void *) indices, sizeof(int), nent, fp) ;
   break ;
case SUBMTX_SPARSE_TRIPLES :
   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   rc = fwrite((void *) rowids, sizeof(int), nent, fp) ;
   rc = fwrite((void *) colids, sizeof(int), nent, fp) ;
   break ;
case SUBMTX_DENSE_SUBROWS :
   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   rc = fwrite((void *) firstlocs, sizeof(int), nrow, fp) ;
   rc = fwrite((void *) sizes,  sizeof(int), nrow, fp) ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                         &firstlocs, &sizes, &entries) ;
   rc = fwrite((void *) firstlocs, sizeof(int), ncol, fp) ;
   rc = fwrite((void *) sizes,  sizeof(int), ncol, fp) ;
   break ;
case SUBMTX_DIAGONAL :
   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
case SUBMTX_BLOCK_DIAGONAL_HERM :
   SubMtx_blockDiagonalInfo(mtx, &nrow, &nent, &pivotsizes, &entries) ;
   rc = fwrite((void *) pivotsizes, sizeof(int), nrow, fp) ;
   break ;
}
/*
   ---------------------
   write out the entries
   ---------------------
*/
if ( SUBMTX_IS_REAL(mtx) ) {
   rc = fwrite((void *) entries, sizeof(double), nent, fp) ;
} else if ( SUBMTX_IS_COMPLEX(mtx) ) {
   rc = fwrite((void *) entries, sizeof(double), 2*nent, fp) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- write a SubMtx object in human readable form 
  
   created -- 98feb07, cca
   -------------------------------------------------------
*/
int
SubMtx_writeForHumanEye (
   SubMtx   *mtx,
   FILE     *fp
) {
A2       a2 ;
double   imag, real ;
int      ierr, irow, jcol, ncol, nrow ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeForHumanEye(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
SubMtx_writeStats(mtx, fp) ;
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
fprintf(fp, "\n rowids : ") ;
IVfp80(fp, nrow, rowind, 80, &ierr) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
fprintf(fp, "\n colids : ") ;
IVfp80(fp, ncol, colind, 80, &ierr) ;

A2_setDefaultFields(&a2) ;
A2_init(&a2, mtx->type, nrow, ncol, 1, nrow, NULL) ;
A2_zero(&a2) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      if ( SUBMTX_IS_REAL(mtx) ) {
         SubMtx_realEntry(mtx, irow, jcol, &real) ;
         A2_setRealEntry(&a2, irow, jcol, real) ;
      } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
         SubMtx_complexEntry(mtx, irow, jcol, &real, &imag) ;
         A2_setComplexEntry(&a2, irow, jcol, real, imag) ;
      }
   }
}
A2_writeForHumanEye(&a2, fp) ;
A2_clearData(&a2) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- write the header and scalar quantities for a SubMtx object
  
   created -- 98feb07, cca
   -------------------------------------------------------------------
*/
int
SubMtx_writeStats (
   SubMtx   *mtx,
   FILE   *fp
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeStats(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
fprintf(fp, 
        "\n\n SubMtx object : type %d, (rowid,colid) = (%d,%d)"
        "\n             : %d rows x %d columns, %d entries"
        "\n             : %d bytes in workspace, %d used, base %p",
        mtx->type, mtx->rowid, mtx->colid,
        mtx->nrow, mtx->ncol, mtx->nent,
        SubMtx_nbytesInWorkspace(mtx), SubMtx_nbytesInUse(mtx),
        SubMtx_workspace(mtx)) ;
switch ( mtx->type ) {
case SPOOLES_REAL :
   fprintf(fp, "\n             : real entries") ;
   break ;
case SPOOLES_COMPLEX :
   fprintf(fp, "\n             : complex entries") ;
   break ;
default :
   fprintf(fp, "\n             : unknown entries") ;
   break ;
}
switch ( mtx->mode ) {
case SUBMTX_DENSE_ROWS :
   fprintf(fp, "\n             : dense storage via rows") ;
   break ;
case SUBMTX_DENSE_COLUMNS :
   fprintf(fp, "\n             : dense storage via columns") ;
   break ;
case SUBMTX_SPARSE_ROWS :
   fprintf(fp, "\n             : sparse storage via rows") ;
   break ;
case SUBMTX_SPARSE_COLUMNS :
   fprintf(fp, "\n             : sparse storage via columns") ;
   break ;
case SUBMTX_SPARSE_TRIPLES :
   fprintf(fp, "\n             : sparse storage via triples") ;
   break ;
case SUBMTX_DENSE_SUBROWS :
   fprintf(fp, "\n             : sparse storage via dense subrows") ;
   break ;
case SUBMTX_DENSE_SUBCOLUMNS :
   fprintf(fp, "\n             : sparse storage via dense subcolumns") ;
   break ;
case SUBMTX_DIAGONAL :
   fprintf(fp, "\n             : diagonal matrix") ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
   fprintf(fp, "\n             : block diagonal symmetric matrix") ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_HERM :
   fprintf(fp, "\n             : block diagonal hermitian matrix") ;
   break ;
default :
   fprintf(fp, "\n             : unknown storage mode") ;
   break ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- write the matrix entries out for matlab
  
   created -- 98feb07, cca
   --------------------------------------------------
*/
void
SubMtx_writeForMatlab (
   SubMtx   *mtx,
   char     *mtxname,
   FILE     *fp
) {
int   ncol, nrow ;
int   *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || mtxname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", mtx, mtxname, fp) ;
   exit(-1) ;
}
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
if ( SUBMTX_IS_DENSE_ROWS(mtx) || SUBMTX_IS_DENSE_COLUMNS(mtx) ) {
   int      ij, inc1, inc2, irow, jcol, ncol, nrow ;
   double   *entries ;

   SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &entries) ;
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
         ij = irow*inc1 + jcol*inc2 ;
         if ( SUBMTX_IS_REAL(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e ; ", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1, entries[ij]) ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1,
                    entries[2*ij], entries[2*ij+1]) ;
         }
      }
   }
} else if ( SUBMTX_IS_SPARSE_ROWS(mtx) ) {
   double   *entries ;
   int      ii, irow, jcol, nent, nrow, rowsize, *indices, *sizes ;

   SubMtx_sparseRowsInfo(mtx, 
                         &nrow, &nent, &sizes, &indices, &entries) ;
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      rowsize = sizes[irow] ;
      for ( ii = 0 ; ii < rowsize ; ii++ ) {
         jcol = indices[ii] ;
         if ( SUBMTX_IS_REAL(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1, entries[ii]) ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1, 
                    entries[2*ii], entries[2*ii+1]) ;
         }
      }
      indices += rowsize ;
      if ( SUBMTX_IS_REAL(mtx) ) {
         entries += rowsize ;
      } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
         entries += 2*rowsize ;
      }
   }
} else if ( SUBMTX_IS_SPARSE_COLUMNS(mtx) ) {
   double   *entries ;
   int      colsize, ii, irow, jcol, nent, ncol, *indices, *sizes ;

   SubMtx_sparseColumnsInfo(mtx, &ncol, &nent, 
                          &sizes, &indices, &entries) ;
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      colsize = sizes[jcol] ;
      for ( ii = 0 ; ii < colsize ; ii++ ) {
         irow = indices[ii] ;
         if ( SUBMTX_IS_REAL(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1, entries[ii]) ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                    rowind[irow] + 1, colind[jcol] + 1, 
                    entries[2*ii], entries[2*ii+1]) ;
         }
      }
      indices += colsize ;
      if ( SUBMTX_IS_REAL(mtx) ) {
         entries += colsize ;
      } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
         entries += 2*colsize ;
      }
   }
} else if ( SUBMTX_IS_SPARSE_TRIPLES(mtx) ) {
   double   *entries ;
   int      ii, irow, jcol, nent, *colids, *rowids ;

   SubMtx_sparseTriplesInfo(mtx, &nent, &rowids, &colids, &entries) ;
   for ( ii = 0 ; ii < nent ; ii++ ) {
      irow = rowids[ii] ;
      jcol = colids[ii] ;
      if ( SUBMTX_IS_REAL(mtx) ) {
         fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname,
                 rowind[irow] + 1, colind[jcol] + 1, entries[ii]) ;
      } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
         fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                 rowind[irow] + 1, colind[jcol] + 1, 
                 entries[2*ii], entries[2*ii+1]) ;
      }
   }
} else if ( SUBMTX_IS_DENSE_SUBROWS(mtx) ) {
   double   *entries ;
   int      first, ii, irow, jcol, last, nent, nrow ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubrowsInfo(mtx, &nrow, &nent, 
                         &firstlocs, &sizes, &entries) ;
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      if ( sizes[irow] > 0 ) {
         first = firstlocs[irow] ;
         last  = first + sizes[irow] - 1 ;
         for ( jcol = first, ii = 0 ; jcol <= last ; jcol++, ii++ ) {
            if ( SUBMTX_IS_REAL(mtx) ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e ;", 
                       mtxname, rowind[irow] + 1, colind[jcol] + 1, 
                       entries[ii]) ;
            } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", 
                       mtxname, rowind[irow] + 1, colind[jcol] + 1, 
                       entries[2*ii], entries[2*ii+1]) ;
            } 
         }
         if ( SUBMTX_IS_REAL(mtx) ) {
            entries += sizes[irow] ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            entries += 2*sizes[irow] ;
         }
      }
   }
} else if ( SUBMTX_IS_DENSE_SUBCOLUMNS(mtx) ) {
   double   *entries ;
   int      first, ii, irow, jcol, last, ncol, nent ;
   int      *firstlocs, *sizes ;

   SubMtx_denseSubcolumnsInfo(mtx, &ncol, &nent, 
                            &firstlocs, &sizes, &entries) ;
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      if ( sizes[jcol] > 0 ) {
         first = firstlocs[jcol] ;
         last  = first + sizes[jcol] - 1 ;
         for ( irow = first, ii = 0 ; irow <= last ; irow++, ii++ ) {
            if ( SUBMTX_IS_REAL(mtx) ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname,
                       rowind[irow] + 1, colind[jcol] + 1, entries[ii]);
            } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
               fprintf(fp, 
                       "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                       rowind[irow] + 1, colind[jcol] + 1, 
                       entries[2*ii], entries[2*ii+1]) ;
            }
         }
         if ( SUBMTX_IS_REAL(mtx) ) {
            entries += sizes[jcol] ;
         } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
            entries += 2*sizes[jcol] ;
         }
      }
   }
} else if ( SUBMTX_IS_DIAGONAL(mtx) ) {
   double   *entries ;
   int      irow, nent ;

   SubMtx_diagonalInfo(mtx, &nent, &entries) ;
   for ( irow = 0 ; irow < nent ; irow++ ) {
      if ( SUBMTX_IS_REAL(mtx) ) {
         fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname,
                 rowind[irow]+1, colind[irow]+1, entries[irow]) ;
      } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
         fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", mtxname,
                 rowind[irow]+1, colind[irow]+1, 
                 entries[2*irow], entries[2*irow+1]) ;
      }
   }
} else if ( SUBMTX_IS_BLOCK_DIAGONAL_SYM(mtx) ) {
   double   *entries ;
   int      ii, ipivot, irow, jj, m, kk, ncol, nent ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   for ( irow = ipivot = kk = 0 ; irow < ncol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++ ) {
         for ( jj = ii ; jj < m ; jj++, kk++ ) {
            if ( SUBMTX_IS_REAL(mtx) ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname, 
                       rowind[irow+ii]+1, colind[irow+jj]+1, 
                       entries[kk]) ;
            } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", 
                       mtxname, rowind[irow + ii]+1, colind[irow+jj]+1, 
                       entries[2*kk], entries[2*kk+1]) ;
            }
            if ( ii != jj ) {
               if ( SUBMTX_IS_REAL(mtx) ) {
                  fprintf(fp, "\n %s(%d,%d) = %20.12e ;", mtxname, 
                          colind[irow + jj]+1, rowind[irow + ii]+1, 
                          entries[kk]) ;
               } else if ( SUBMTX_IS_COMPLEX(mtx) ) {
                  fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", 
                          mtxname, colind[irow + jj]+1, 
                          rowind[irow + ii]+1, 
                          entries[2*kk], entries[2*kk+1]) ;
               }
            }
         }
      }
      irow += m ;
   }
} else if ( SUBMTX_IS_BLOCK_DIAGONAL_HERM(mtx) ) {
   double   *entries ;
   int      ii, ipivot, irow, jj, m, kk, ncol, nent ;
   int      *pivotsizes ;

   SubMtx_blockDiagonalInfo(mtx, &ncol, &nent, &pivotsizes, &entries) ;
   for ( irow = ipivot = kk = 0 ; irow < ncol ; ipivot++ ) {
      m = pivotsizes[ipivot] ;
      for ( ii = 0 ; ii < m ; ii++ ) {
         for ( jj = ii ; jj < m ; jj++, kk++ ) {
            fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", 
                    mtxname, rowind[irow+ii]+1, colind[irow+jj]+1, 
                    entries[2*kk], entries[2*kk+1]) ;
            if ( ii != jj ) {
               fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;", 
                       mtxname, colind[irow+jj]+1, rowind[irow+ii]+1,
                       entries[2*kk], -entries[2*kk+1]) ;
            }
         }
      }
      irow += m ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
