/*  IO.c  */

#include "../DenseMtx.h"

static const char *suffixb = ".densemtxb" ;
static const char *suffixf = ".densemtxf" ;

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to read in the object from a file

   input --

      fn -- filename, must be *.densemtxb or *.densemtxf

   return value -- 1 if success, 0 if failure
 
   created -- 98aug13
   -----------------------------------------------------------
*/
int
DenseMtx_readFromFile ( 
   DenseMtx   *mtx, 
   char       *fn 
) {
FILE   *fp ;
int    fnlength, rc = 0, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_readFromFile(%p,%s)"
           "\n bad input", mtx, fn) ;
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
         fprintf(stderr, "\n error in DenseMtx_readFromFile()"
                 "\n unable to open file %s", fn) ;
      } else {
         rc = DenseMtx_readFromBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in DenseMtx_readFromFile()"
                 "\n unable to open file %s", fn) ;
      } else {
         rc = DenseMtx_readFromFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in DenseMtx_readFromFile()"
              "\n bad DenseMtx file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, 
           "\n error in DenseMtx_readFromFile()"
           "\n bad DenseMtx file name %s,"
           "\n must end in %s (binary) or %s (formatted)\n",
           fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to read an object from a formatted file

   return value -- 1 if success, 0 if failure
 
   created -- 98aug13
   --------------------------------------------------
*/
int
DenseMtx_readFromFormattedFile ( 
   DenseMtx   *mtx, 
   FILE       *fp 
) {
int   rc, size ;
int   itemp[7] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile(%p,%p)"
           "\n bad input", mtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
DenseMtx_clearData(mtx) ;
/*
   --------------------------------------------
   read in the seven scalar parameters: 
     type, rowid, colid, nrow, ncol, inc1, inc2
   --------------------------------------------
*/
if ( (rc = IVfscanf(fp, 7, itemp)) != 7 ) {
   fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile()"
           "\n %d items of %d read\n", rc, 7) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
DenseMtx_init(mtx, itemp[0], itemp[1], itemp[2], 
              itemp[3], itemp[4], itemp[5], itemp[6]) ;
/*
   ---------------------------
   read in the rowind[] vector
   ---------------------------
*/
if ( (size = mtx->nrow) > 0 ) {
   if ( (rc = IVfscanf(fp, size, mtx->rowind)) != size ) {
      fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile()"
              "\n %d items of %d read for rowind\n", rc, size) ;
      return(0) ;
   }
}
/*
   ---------------------------
   read in the colind[] vector
   ---------------------------
*/
if ( (size = mtx->ncol) > 0 ) {
   if ( (rc = IVfscanf(fp, size, mtx->colind)) != size ) {
      fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile()"
              "\n %d items of %d read for colind\n", rc, size) ;
      return(0) ;
   }
}
/*
   ----------------------------
   read in the entries[] vector
   ----------------------------
*/
if ( (size = mtx->nrow*mtx->ncol) > 0 ) {
   if ( DENSEMTX_IS_REAL(mtx) ) {
      if ( (rc = DVfscanf(fp, size, mtx->entries)) != size ) {
         fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile()"
                 "\n %d items of %d read for entries\n", rc, size) ;
         return(0) ;
      }
   } else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
      if ( (rc = DVfscanf(fp, 2*size, mtx->entries)) != 2*size ) {
         fprintf(stderr, "\n error in DenseMtx_readFromFormattedFile()"
                 "\n %d items of %d read for entries\n", rc, 2*size) ;
         return(0) ;
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to read an object from a binary file

   return value -- 1 if success, 0  if failure
 
   created -- 98aug13
   -----------------------------------------------
*/
int
DenseMtx_readFromBinaryFile ( 
   DenseMtx   *mtx, 
   FILE       *fp 
) {
int   rc, size ;
int   itemp[7] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DenseMtx_readFromBinaryFile(%p,%p)"
           "\n bad input", mtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
DenseMtx_clearData(mtx) ;
/*
   --------------------------------------------
   read in the seven scalar parameters: 
     type, rowid, colid, nrow, ncol, inc1, inc2
   --------------------------------------------
*/
if ( (rc = fread((char *) itemp, sizeof(int), 7, fp)) != 7 ) {
   fprintf(stderr, "\n error in DenseMtx_readFromBinaryFile()"
           "\n %d items of %d read\n", rc, 7) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
DenseMtx_init(mtx, itemp[0], itemp[1], itemp[2], 
              itemp[3], itemp[4], itemp[5], itemp[6]) ;
/*
   ---------------------------
   read in the rowind[] vector
   ---------------------------
*/
if ( (size = mtx->nrow) > 0 ) {
   if ( (rc = fread(mtx->rowind, sizeof(int), size, fp)) != size ) {
      fprintf(stderr, "\n error in DenseMtx_readFromBinaryFile()"
              "\n %d items of %d read for rowind[]\n", rc, size) ;
      return(0) ;
   }
}
/*
   ---------------------------
   read in the colind[] vector
   ---------------------------
*/
if ( (size = mtx->ncol) > 0 ) {
   if ( (rc = fread(mtx->colind, sizeof(int), size, fp)) != size ) {
      fprintf(stderr, "\n error in DenseMtx_readFromBinaryFile()"
              "\n %d items of %d read for colind[]\n", rc, size) ;
      return(0) ;
   }
}
/*
   ----------------------------
   read in the entries[] vector
   ----------------------------
*/
if ( (size = mtx->nrow*mtx->ncol) > 0 ) {
   if ( DENSEMTX_IS_REAL(mtx) ) {
      rc = fread(mtx->entries, sizeof(double), size, fp) ;
      if ( rc != size ) {
         fprintf(stderr, "\n error in DenseMtx_readFromBinaryFile()"
                 "\n %d items of %d read for entries\n", rc, size) ;
         return(0) ;
      }
   } else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
      rc = fread(mtx->entries, sizeof(double), 2*size, fp) ;
      if ( rc != 2*size ) {
         fprintf(stderr, "\n error in DenseMtx_readFromBinaryFile()"
                 "\n %d items of %d read for entries\n", rc, 2*size) ;
         return(0) ;
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   purpose -- to write an object to a file

   input --

      fn -- filename
        *.densemtxb -- binary
        *.densemtxf -- formatted
        anything else -- for human eye

   error return --
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
     -3 -- unable to open file
 
   created -- 98aug13
   ---------------------------------------
*/
int
DenseMtx_writeToFile ( 
   DenseMtx   *mtx, 
   char       *fn 
) {
FILE   *fp ;
int    fnlength, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToFile(%p,%s)"
           "\n mtx is NULL", mtx, fn) ;
   return(-1) ;
}
if ( fn == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToFile(%p,%s)"
           "\n fn is NULL", mtx, fn) ;
   return(-2) ;
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
         fprintf(stderr, "\n error in DenseMtx_writeToFile()"
                 "\n unable to open file %s", fn) ;
         return(-3) ;
      } else {
         DenseMtx_writeToBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in DenseMtx_writeToFile()"
                 "\n unable to open file %s", fn) ;
         return(-3) ;
      } else {
         DenseMtx_writeToFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in DenseMtx_writeToFile()"
                 "\n unable to open file %s", fn) ;
         return(-3) ;
      } else {
         DenseMtx_writeForHumanEye(mtx, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in DenseMtx_writeToFile()"
              "\n unable to open file %s", fn) ;
      return(-3) ;
   } else {
      DenseMtx_writeForHumanEye(mtx, fp) ;
      fclose(fp) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an object to a formatted file

   return value -- 
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
 
   created -- 98aug13
   -------------------------------------------------
*/
int
DenseMtx_writeToFormattedFile ( 
   DenseMtx   *mtx, 
   FILE       *fp 
) {
int   size ;

if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToFormattedFile()"
           "\n mtx is NULL") ;
   return(-1) ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToFormattedFile()"
           "\n fp is NULL") ;
   return(-2) ;
}
fprintf(fp, "\n %d %d %d %d %d %d %d", 
        mtx->type, mtx->rowid, mtx->colid, 
        mtx->nrow, mtx->ncol, mtx->inc1, mtx->inc2) ;
if ( (size = mtx->nrow) > 0 && mtx->rowind != NULL ) {
   IVfprintf(fp, size, mtx->rowind) ;
}
if ( (size = mtx->ncol) > 0 && mtx->colind != NULL ) {
   IVfprintf(fp, size, mtx->colind) ;
}
if ( (size = mtx->nrow*mtx->ncol) > 0 && mtx->entries != NULL ) {
   if ( DENSEMTX_IS_REAL(mtx) ) {
      DVfprintf(fp, size, mtx->entries) ;
   } else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
      DVfprintf(fp, 2*size, mtx->entries) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to write an adjacency object to a binary file

   return value -- 
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fn is NULL
 
   created -- 98aug13
   --------------------------------------------------------
*/
int
DenseMtx_writeToBinaryFile ( 
   DenseMtx   *mtx, 
   FILE       *fp 
) {
int   size ;

if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToBinaryFile()"
           "\n mtx is NULL") ;
   return(-1) ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeToBinaryFile()"
           "\n fp is NULL") ;
   return(-2) ;
}
fwrite((void *) &mtx->type,  sizeof(int), 1, fp) ;
fwrite((void *) &mtx->rowid, sizeof(int), 1, fp) ;
fwrite((void *) &mtx->colid, sizeof(int), 1, fp) ;
fwrite((void *) &mtx->nrow,  sizeof(int), 1, fp) ;
fwrite((void *) &mtx->ncol,  sizeof(int), 1, fp) ;
fwrite((void *) &mtx->inc1,  sizeof(int), 1, fp) ;
fwrite((void *) &mtx->inc2,  sizeof(int), 1, fp) ;
if ( (size = mtx->nrow) > 0 && mtx->rowind != NULL ) {
   fwrite((void *) mtx->rowind, sizeof(int), size, fp) ;
}
if ( (size = mtx->ncol) > 0 && mtx->colind != NULL ) {
   fwrite((void *) mtx->colind, sizeof(int), size, fp) ;
}
if ( (size = mtx->nrow*mtx->ncol) > 0 && mtx->entries != NULL ) {
   if ( DENSEMTX_IS_REAL(mtx) ) {
      fwrite((void *) mtx->entries, sizeof(double), size, fp) ;
   } else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
      fwrite((void *) mtx->entries, sizeof(double), 2*size, fp) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write the object's statistics to a file
              in human readable form

   return value -- 
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fp is NULL

   created -- 98may02, cca
   -----------------------------------------------------
*/
int
DenseMtx_writeStats (
   DenseMtx   *mtx,
   FILE       *fp
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeStats()"
           "\n mtx is NULL") ;
   return(-1) ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeStats()"
           "\n fp is NULL") ;
   return(-2) ;
}
fprintf(fp, "\n DenseMtx object at address %p", mtx) ;
switch ( mtx->type ) {
case SPOOLES_REAL :
   fprintf(fp, ", real entries") ;
   break ;
case SPOOLES_COMPLEX :
   fprintf(fp, ", complex entries") ;
   break ;
default :
   fprintf(fp, ", unknown entries type") ;
   break ;
}
fprintf(fp, "\n row id = %d, col id = %d"
        "\n nrow = %d, ncol = %d, inc1 = %d, inc2 = %d",
        mtx->rowid, mtx->colid, 
        mtx->nrow, mtx->ncol, mtx->inc1, mtx->inc2) ;
fprintf(fp, "\n rowind = %p, colind = %p, entries = %p",
        mtx->rowind, mtx->colind, mtx->entries) ;
fprintf(fp, ", base = %p", DV_entries(&mtx->wrkDV)) ;
fprintf(fp, 
       "\n rowind - base = %d, colind - base = %d, entries - base = %d",
       mtx->rowind - (int *) DV_entries(&mtx->wrkDV),
       mtx->colind - (int *) DV_entries(&mtx->wrkDV),
       mtx->entries - DV_entries(&mtx->wrkDV)) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   return value -- 
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- fp is NULL

   created -- 98may02, cca
   ----------------------------------------
*/
int
DenseMtx_writeForHumanEye (
   DenseMtx   *mtx,
   FILE       *fp
) {
A2    a2 ;
int   ierr, ncol, nrow ;
int   *colind, *rowind ; 
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeForHumanEye()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeForHumanEye()"
           "\n mtx is NULL\n") ;
   return(-2) ;
}
DenseMtx_writeStats(mtx, fp) ;
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
if ( nrow > 0 && rowind != NULL ) {
   fprintf(fp, "\n mtx's row indices at %p", rowind) ;
   IVfp80(fp, nrow, rowind, 80, &ierr) ;
}
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
if ( ncol > 0 && colind != NULL ) {
   fprintf(fp, "\n mtx's column indices at %p", colind) ;
   IVfp80(fp, ncol, colind, 80, &ierr) ;
}
if ( nrow > 0 && ncol > 0 ) {
   A2_setDefaultFields(&a2) ;
   DenseMtx_setA2(mtx, &a2) ;
   A2_writeForHumanEye(&a2, fp) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to write the object to a matlab file

   return value -- 
      1 -- normal return
     -1 -- mtx is NULL
     -2 -- mtx is NULL
     -3 -- fp is NULL

   created -- 98may02, cca
   -----------------------------------------------
*/
int
DenseMtx_writeForMatlab (
   DenseMtx   *mtx,
   char       *mtxname,
   FILE       *fp
) {
double   *entries ;
int      inc1, inc2, irow, jcol, ncol, nrow ;
int      *colind, *rowind ; 
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeForMatlab()"
           "\n mtx is NULL\n") ;
   return(-1) ;
}
if ( mtxname == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeForMatlab()"
           "\n mtxname is NULL\n") ;
   return(-2) ;
}
if ( fp == NULL ) {
   fprintf(stderr, "\n fatal error in DenseMtx_writeForMatlab()"
           "\n fp is NULL\n") ;
   return(-3) ;
}
DenseMtx_rowIndices(mtx, &nrow, &rowind) ;
DenseMtx_columnIndices(mtx, &ncol, &colind) ;
DenseMtx_dimensions(mtx, &nrow, &ncol) ;
inc1 = DenseMtx_rowIncrement(mtx) ;
inc2 = DenseMtx_columnIncrement(mtx) ;
entries = DenseMtx_entries(mtx) ;
if ( DENSEMTX_IS_REAL(mtx) ) {
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                 mtxname, rowind[irow]+1, colind[jcol]+1,
                 entries[irow*inc1+jcol*inc2]) ;
      }
   }
} else if ( DENSEMTX_IS_COMPLEX(mtx) ) {
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      for ( irow = 0 ; irow < nrow ; irow++ ) {
         fprintf(fp, "\n %s(%d,%d) = %24.16e + %24.16e*i ;",
                 mtxname, rowind[irow]+1, colind[jcol]+1,
                 entries[2*(irow*inc1+jcol*inc2)],
                 entries[2*(irow*inc1+jcol*inc2)+1]) ;
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
