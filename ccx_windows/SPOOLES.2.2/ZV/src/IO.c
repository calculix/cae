/*  IO.c  */

#include "../ZV.h"

static const char *suffixb = ".zvb" ;
static const char *suffixf = ".zvf" ;

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to read in an ZV object from a file

   input --

      fn -- filename, must be *.zvb or *.zvf

   return value -- 1 if success, 0 if failure

   created -- 98jan22, cca
   ----------------------------------------------
*/
int
ZV_readFromFile ( 
   ZV     *zv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in ZV_readFromFile(%p,%s)"
    "\n bad input\n", zv, fn) ;
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
         fprintf(stderr, "\n error in ZV_readFromFile(%p,%s)"
                 "\n unable to open file %s", zv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ZV_readFromBinaryFile(zv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in ZV_readFromFile(%p,%s)"
                 "\n unable to open file %s", zv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ZV_readFromFormattedFile(zv, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in ZV_readFromFile(%p,%s)"
              "\n bad ZV file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              zv, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in ZV_readFromFile(%p,%s)"
       "\n bad ZV file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       zv, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to read an ZV object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 98jan22, cca
   -----------------------------------------------------
*/
int
ZV_readFromFormattedFile ( 
   ZV     *zv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in ZV_readFromFormattedFile(%p,%p)"
           "\n bad input\n", zv, fp) ;
   return(0) ;
}
ZV_clearData(zv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fscanf(fp, "%d", &size)) != 1 ) {
   fprintf(stderr, "\n error in ZV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", zv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
ZV_init(zv, size, NULL) ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = DVfscanf(fp, 2*size, ZV_entries(zv))) != 2*size ) {
   fprintf(stderr, "\n error in ZV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", zv, fp, rc, 2*size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read an ZV object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 98jan22, cca
   ---------------------------------------------------
*/
int
ZV_readFromBinaryFile ( 
   ZV    *zv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_readFromBinaryFile(%p,%p)"
           "\n bad input\n", zv, fp) ;
   return(0) ;
}
ZV_clearData(zv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fread((void *) &size, sizeof(int), 1, fp)) != 1 ) {
   fprintf(stderr, "\n error in ZV_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", zv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
ZV_init(zv, size, NULL) ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = fread((void *) ZV_entries(zv), sizeof(double), 2*size, fp)) 
      != 2*size ) {
   fprintf(stderr, "\n error in ZV_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", 
           zv, fp, rc, 2*size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write an ZV object to a file

   input --

      fn -- filename
        *.zvb -- binary
        *.zvf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 98jan22, cca
   -------------------------------------------
*/
int
ZV_writeToFile ( 
   ZV    *zv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_writeToFile(%p,%s)"
    "\n bad input\n", zv, fn) ; 
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
         fprintf(stderr, "\n error in ZV_writeToFile(%p,%s)"
                 "\n unable to open file %s", zv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ZV_writeToBinaryFile(zv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in ZV_writeToFile(%p,%s)"
                 "\n unable to open file %s", zv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ZV_writeToFormattedFile(zv, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in ZV_writeToFile(%p,%s)"
                 "\n unable to open file %s", zv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ZV_writeForHumanEye(zv, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in ZV_writeToFile(%p,%s)"
              "\n unable to open file %s", zv, fn, fn) ;
      rc = 0 ;
   } else {
      rc = ZV_writeForHumanEye(zv, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an ZV object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 98jan22, cca
   -----------------------------------------------------
*/
int
ZV_writeToFormattedFile ( 
   ZV    *zv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fp == NULL || zv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in ZV_writeToFormattedFile(%p,%p)"
           "\n bad input\n", zv, fp) ;
   fprintf(stderr, "\n zv->size = %d", zv->size) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the size of the vector
   -------------------------------------
*/
size = ZV_size(zv) ;
rc   = fprintf(fp, "\n %d", size) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ZV_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", zv, fp, rc) ;
   return(0) ;
}
if ( size > 0 ) {
   DVfprintf(fp, 2*size, ZV_entries(zv)) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an ZV object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 98jan22, cca
   --------------------------------------------------
*/
int
ZV_writeToBinaryFile ( 
   ZV    *zv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fp == NULL || zv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in ZV_writeToBinaryFile(%p,%p)"
           "\n bad input\n", zv, fp) ;
   exit(-1) ;
}
size = ZV_size(zv) ;
rc = fwrite((void *) &size, sizeof(int), 1, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in ZV_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", zv, fp, rc, 1) ;
   return(0) ;
}
rc = fwrite((void *) ZV_entries(zv), sizeof(double), 2*size, fp) ;
if ( rc != 2*size ) {
   fprintf(stderr, "\n error in ZV_writeToBinaryFile(%p,%p)"
           "\n %d of %d items written\n",
           zv, fp, rc, 2*size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an ZV object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 98jan22, cca
   -------------------------------------------------
*/
int
ZV_writeForHumanEye ( 
   ZV    *zv, 
   FILE   *fp 
) {
double   *vec ;
int      ii, jj, rc, size ;

if ( zv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_writeForHumanEye(%p,%p)"
           "\n bad input\n", zv, fp) ;
   exit(-1) ;
}
if ( (rc = ZV_writeStats(zv, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in ZV_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from ZV_writeStats(%p,%p)\n",
           zv, fp, rc, zv, fp) ;
   return(0) ;
}
size = ZV_size(zv) ;
vec  = ZV_entries(zv) ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   if ( ii % 2 == 0 ) {
      fprintf(fp, "\n") ;
   }
   fprintf(fp, " < %12.4e, %12.4e >", vec[jj], vec[jj+1]) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the ZV object

   return value -- 1 if success, 0 otherwise

   created -- 98jan22, cca
   ---------------------------------------------------------
*/
int
ZV_writeStats ( 
   ZV    *zv, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in ZV_writeStats(%p,%p)"
           "\n bad input\n", zv, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n ZV : double complex vector object : ") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, 
             " size = %d, maxsize = %d, owned = %d", 
             zv->size, zv->maxsize, zv->owned) ;
if ( rc < 0 ) { goto IO_error ; }
return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in ZV_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", zv, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- write the vector entries out for matlab
  
   created -- 98apr15, cca
   --------------------------------------------------
*/
void
ZV_writeForMatlab (
   ZV     *zv,
   char   *vecname,
   FILE   *fp
) {
int      ii, jj, size ;
double   *z ;
/*
   ---------------
   check the input
   ---------------
*/
if ( zv == NULL || vecname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in ZV_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", zv, vecname, fp) ;
   exit(-1) ;
}
ZV_sizeAndEntries(zv, &size, &z) ;
for ( ii = jj = 0 ; ii < size ; ii++, jj += 2 ) {
   fprintf(fp, "\n %s(%d) = %24.16e + %24.16e*i;", 
           vecname, ii + 1, z[jj], z[jj+1]) ;
}
return ; }

/*--------------------------------------------------------------------*/
