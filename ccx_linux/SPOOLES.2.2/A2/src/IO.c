/*  IO.c  */

#include "../A2.h"

static const char *suffixb = ".a2b" ;
static const char *suffixf = ".a2f" ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to read in the object from a file

   input --

      fn -- filename, must be *.a2b or *.a2f

   return value -- 1 if success, 0 if failure

   created -- 98may01, cca
   --------------------------------------------
*/
int
A2_readFromFile ( 
   A2     *mtx, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc = 0, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in A2_readFromFile(%p,%s)"
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
         fprintf(stderr, "\n error in A2_readFromFile(%s)"
                 "\n unable to open file %s", fn, fn) ;
      } else {
         rc = A2_readFromBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in A2_readFromFile(%s)"
                 "\n unable to open file %s", fn, fn) ;
      } else {
         rc = A2_readFromFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in A2_readFromFile(%s)"
              "\n bad A2 file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in A2_readFromFile(%s)"
           "\n bad A2 file name %s,"
           "\n must end in %s (binary) or %s (formatted)\n",
           fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to read an object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 98may01, cca
   --------------------------------------------------
*/
int
A2_readFromFormattedFile ( 
   A2    *mtx, 
   FILE   *fp 
) {
int   rc, size ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in A2_readFromFormattedFile(%p,%p)"
           "\n bad input", mtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
A2_clearData(mtx) ;
/*
   -----------------------------------------------------------
   read in the five scalar parameters: type n1, n2, inc1, inc2
   -----------------------------------------------------------
*/
if ( (rc = IVfscanf(fp, 5, itemp)) != 5 ) {
   fprintf(stderr, "\n error in A2_readFromFormattedFile()"
           "\n %d items of %d read\n", rc, 5) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
A2_init(mtx, itemp[0], itemp[1], itemp[2], itemp[3], itemp[4], NULL) ;
/*
   ----------------------------
   read in the entries[] vector
   ----------------------------
*/
if ( (size = 1 + (mtx->n1-1)*mtx->inc1 + (mtx->n2-1)*mtx->inc2) > 0 ) {
   if ( A2_IS_REAL(mtx) ) {
      if ( (rc = DVfscanf(fp, size, mtx->entries)) != size ) {
         fprintf(stderr, "\n error in A2_readFromFormattedFile"
                 "\n %d items of %d read\n", rc, size) ;
         return(0) ;
      }
   } else if ( A2_IS_COMPLEX(mtx) ) {
      if ( (rc = DVfscanf(fp, 2*size, mtx->entries)) != 2*size ) {
         fprintf(stderr, "\n error in A2_readFromFormattedFile"
                 "\n %d items of %d read\n", rc, 2*size) ;
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

   created -- 98may01, cca
   -----------------------------------------------
*/
int
A2_readFromBinaryFile ( 
   A2    *mtx, 
   FILE   *fp 
) {
int   rc, size ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in A2_readFromBinaryFile(%p,%p)"
           "\n bad input", mtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
A2_clearData(mtx) ;
/*
   ------------------------------------------------------------
   read in the five scalar parameters: type, n1, n2, inc1, inc2
   ------------------------------------------------------------
*/
if ( (rc = fread((char *) itemp, sizeof(int), 5, fp)) != 5 ) {
   fprintf(stderr, "\n error in A2_readFromBinaryFile"
           "\n %d items of %d read\n", rc, 5) ;
   return(0) ;
}
fprintf(stdout, "\n itemp = {%d, %d, %d, %d, %d}", 
        itemp[0], itemp[1], itemp[2], itemp[3], itemp[4]) ;
fflush(stdout) ;
/*
   ---------------------
   initialize the object
   ---------------------
*/
A2_init(mtx, itemp[0], itemp[1], itemp[2], itemp[3], itemp[4], NULL) ;
/*
   ----------------------------
   read in the entries[] vector
   ----------------------------
*/
if ( (size = 1 + (mtx->n1-1)*mtx->inc1 + (mtx->n2-1)*mtx->inc2) > 0 ) {
   if ( A2_IS_REAL(mtx) ) {
      if ( (rc = fread(mtx->entries, sizeof(double), size, fp)) 
           != size ) {
         fprintf(stderr, "\n error in A2_readFromBinaryFile"
                 "\n %d items of %d read\n", rc, size) ;
         return(0) ;
      }
   } else if ( A2_IS_COMPLEX(mtx) ) {
      if ( (rc = fread(mtx->entries, sizeof(double), 2*size, fp)) 
           != 2*size ) {
         fprintf(stderr, "\n error in A2_readFromBinaryFile"
                 "\n %d items of %d read\n", rc, 2*size) ;
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
        *.a2b -- binary
        *.a2f -- formatted
        anything else -- for human eye

   created -- 98may01, cca
   ---------------------------------------
*/
void
A2_writeToFile ( 
   A2    *mtx, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( mtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in A2_writeToFile(%p,%s)"
           "\n bad input", mtx, fn) ;
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
         fprintf(stderr, "\n error in A2_writeToFile()"
                 "\n unable to open file %s", fn) ;
      } else {
         A2_writeToBinaryFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in A2_writeToFile()"
                 "\n unable to open file %s", fn) ;
      } else {
         A2_writeToFormattedFile(mtx, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in A2_writeToFile()"
                 "\n unable to open file %s", fn) ;
      } else {
         A2_writeForHumanEye(mtx, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in A2_writeToFile"
              "\n unable to open file %s", fn) ;
   } else {
      A2_writeForHumanEye(mtx, fp) ;
      fclose(fp) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an object to a formatted file

   created -- 98may01, cca
   -------------------------------------------------
*/
void
A2_writeToFormattedFile ( 
   A2    *mtx, 
   FILE   *fp 
) {
int   size ;

if ( mtx == NULL || fp == NULL ) {
   return ;
}
fprintf(fp, "\n %d %d %d %d %d", 
        mtx->type, mtx->n1, mtx->n2, mtx->inc1, mtx->inc2) ;
if ( (size = 1 + (mtx->n1-1)*mtx->inc1 + (mtx->n2-1)*mtx->inc2) > 0
     && mtx->entries != NULL ) {
   if ( A2_IS_REAL(mtx) ) {
      DVfprintf(fp, size, mtx->entries) ;
   } else if ( A2_IS_COMPLEX(mtx) ) {
      DVfprintf(fp, 2*size, mtx->entries) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to write an adjacency object to a binary file

   created -- 98may01, cca
   --------------------------------------------------------
*/
void
A2_writeToBinaryFile ( 
   A2    *mtx, 
   FILE   *fp 
) {
int   size ;

if ( fp == NULL ) {
   return ;
}
fwrite((void *) &mtx->type, sizeof(int), 1, fp) ;
fwrite((void *) &mtx->n1,   sizeof(int), 1, fp) ;
fwrite((void *) &mtx->n2,   sizeof(int), 1, fp) ;
fwrite((void *) &mtx->inc1, sizeof(int), 1, fp) ;
fwrite((void *) &mtx->inc2, sizeof(int), 1, fp) ;
if ( (size = 1 + (mtx->n1-1)*mtx->inc1 + (mtx->n2-1)*mtx->inc2) > 0
     && mtx->entries != NULL ) {
   if ( A2_IS_REAL(mtx) ) {
      fwrite((void *) &mtx->entries, sizeof(double), size, fp) ;
   } else if ( A2_IS_COMPLEX(mtx) ) {
      fwrite((void *) &mtx->entries, sizeof(double), 2*size, fp) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to write the object for a human eye

   created -- 98may01, cca
   ----------------------------------------------
*/
void
A2_writeForHumanEye ( 
   A2    *mtx, 
   FILE   *fp 
) {
int   i, j, jfirst, jlast, loc ;

if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in A2_writeForHumanEye(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
A2_writeStats(mtx, fp) ;
if ( A2_IS_REAL(mtx) ) {
   for ( jfirst = 0 ; jfirst < mtx->n2 ; jfirst += 4 ) {
      jlast = jfirst + 3 ;
      if ( jlast >= mtx->n2 ) {
         jlast = mtx->n2 - 1 ;
      }
      fprintf(fp, "\n    ") ;
      for ( j = jfirst ; j <= jlast ; j++ ) {
         fprintf(fp, "%19d", j) ;
      }
      for ( i = 0 ; i < mtx->n1 ; i++ ) {
         fprintf(fp, "\n%4d", i) ;
         for ( j = jfirst ; j <= jlast ; j++ ) {
            fprintf(fp, "%19.11e", 
                     mtx->entries[i*mtx->inc1 + j*mtx->inc2]) ;
         }
      }
   }
} else if ( A2_IS_COMPLEX(mtx) ) {
   for ( jfirst = 0 ; jfirst < mtx->n2 ; jfirst += 2 ) {
      jlast = jfirst + 1 ;
      if ( jlast >= mtx->n2 ) {
         jlast = mtx->n2 - 1 ;
      }
      fprintf(fp, "\n    ") ;
      for ( j = jfirst ; j <= jlast ; j++ ) {
         fprintf(fp, "%36d", j) ;
      }
      for ( i = 0 ; i < mtx->n1 ; i++ ) {
         fprintf(fp, "\n%4d", i) ;
         for ( j = jfirst ; j <= jlast ; j++ ) {
            loc = 2*(i*mtx->inc1 + j*mtx->inc2) ;
            fprintf(fp, " (%16.9e,%16.9e*i)", 
                    mtx->entries[loc], mtx->entries[loc+1]) ;
         }
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   purpose -- to write out the statistics

   created -- 98may01, cca
   --------------------------------------
*/
void
A2_writeStats (
   A2    *mtx, 
   FILE   *fp 
) {
if ( mtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in A2_writeStats(%p,%p)"
           "\n bad input\n", mtx, fp) ;
   exit(-1) ;
}
if ( A2_IS_REAL(mtx) ) {
   fprintf(fp, "\n A2 : double 2D array object :") ;
} else if ( A2_IS_COMPLEX(mtx) ) {
   fprintf(fp, "\n A2 : double complex 2D array object :") ;
}
fprintf(fp, 
        "\n %d x %d matrix, inc1 = %d, inc2 = %d,"
        "\n nowned = %d, %d entries allocated, occupies %d bytes" 
        "\n entries = %p, maxabs = %20.12e",
        mtx->n1, mtx->n2, mtx->inc1, mtx->inc2, mtx->nowned,
        mtx->n1*mtx->n2,
        A2_sizeOf(mtx), mtx->entries, A2_maxabs(mtx)) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to write the matrix in matlab format

   created -- 98may01, cca
   -----------------------------------------------
*/
void
A2_writeForMatlab (
   A2    *mtx, 
   char   *mtxname,
   FILE   *fp 
) {
int      irow, jcol, ncol, nrow ;

if ( mtx == NULL || mtxname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in A2_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", mtx, mtxname, fp) ;
   exit(-1) ;
}
nrow = A2_nrow(mtx) ;
ncol = A2_ncol(mtx) ;
for ( irow = 0 ; irow < nrow ; irow++ ) {
   for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
      if ( A2_IS_REAL(mtx) ) {
         double   value ;
         A2_realEntry(mtx, irow, jcol, &value) ;
         fprintf(fp, "\n %s(%d,%d) = %24.16e ;",
                 mtxname, irow+1, jcol+1, value) ;
      } else if ( A2_IS_COMPLEX(mtx) ) {
         double   imag, real ;
         A2_complexEntry(mtx, irow, jcol, &real, &imag) ;
         fprintf(fp, "\n %s(%d,%d) = %24.16e + %24.16e*i ;",
                 mtxname, irow+1, jcol+1, real, imag) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
