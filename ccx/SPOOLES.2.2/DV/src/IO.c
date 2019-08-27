/*  IO.c  */

#include "../DV.h"

static const char *suffixb = ".dvb" ;
static const char *suffixf = ".dvf" ;

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to read in an DV object from a file

   input --

      fn -- filename, must be *.dvb or *.dvf

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   ----------------------------------------------
*/
int
DV_readFromFile ( 
   DV    *dv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in DV_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", dv, fn, __FILE__, __LINE__) ;
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
         fprintf(stderr, "\n error in DV_readFromFile(%p,%s)"
                 "\n unable to open file %s", dv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DV_readFromBinaryFile(dv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in DV_readFromFile(%p,%s)"
                 "\n unable to open file %s", dv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DV_readFromFormattedFile(dv, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in DV_readFromFile(%p,%s)"
              "\n bad DV file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              dv, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in DV_readFromFile(%p,%s)"
       "\n bad DV file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       dv, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to read an DV object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
DV_readFromFormattedFile ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DV_readFromFormattedFile(%p,%p)"
           "\n bad input\n", dv, fp) ;
   return(0) ;
}
DV_clearData(dv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fscanf(fp, "%d", &size)) != 1 ) {
   fprintf(stderr, "\n error in DV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", dv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
DV_init(dv, size, NULL) ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = DVfscanf(fp, size, DV_entries(dv))) != size ) {
   fprintf(stderr, "\n error in DV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", dv, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read an DV object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96jun23, cca
   ---------------------------------------------------
*/
int
DV_readFromBinaryFile ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in DV_readFromBinaryFile(%p,%p)"
           "\n bad input\n", dv, fp) ;
   return(0) ;
}
DV_clearData(dv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fread((void *) &size, sizeof(int), 1, fp)) != 1 ) {
   fprintf(stderr, "\n error in DV_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", dv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
DV_init(dv, size, NULL) ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = fread((void *) DV_entries(dv), sizeof(double), size, fp)) 
      != size ) {
   fprintf(stderr, "\n error in DV_readFromBinaryFile(%p,%p)"
           "\n sizes(%d) : %d items of %d read\n", 
           dv, fp, size, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write an DV object to a file

   input --

      fn -- filename
        *.dvb -- binary
        *.dvf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -------------------------------------------
*/
int
DV_writeToFile ( 
   DV    *dv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in DV_writeToFile(%p,%s)"
    "\n bad input\n", dv, fn) ; 
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
         fprintf(stderr, "\n error in DV_writeToFile(%p,%s)"
                 "\n unable to open file %s", dv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DV_writeToBinaryFile(dv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in DV_writeToFile(%p,%s)"
                 "\n unable to open file %s", dv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DV_writeToFormattedFile(dv, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in DV_writeToFile(%p,%s)"
                 "\n unable to open file %s", dv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DV_writeForHumanEye(dv, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in DV_writeToFile(%p,%s)"
              "\n unable to open file %s", dv, fn, fn) ;
      rc = 0 ;
   } else {
      rc = DV_writeForHumanEye(dv, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an DV object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -----------------------------------------------------
*/
int
DV_writeToFormattedFile ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL || dv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in DV_writeToFormattedFile(%p,%p)"
           "\n bad input\n", dv, fp) ;
   fprintf(stderr, "\n dv->size = %d", dv->size) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the size of the vector
   -------------------------------------
*/
size = DV_size(dv) ;
rc = fprintf(fp, "\n %d", size) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in DV_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", dv, fp, rc) ;
   return(0) ;
}
if ( size > 0 ) {
   DVfprintf(fp, size, DV_entries(dv)) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an DV object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   --------------------------------------------------
*/
int
DV_writeToBinaryFile ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL || dv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in DV_writeToBinaryFile(%p,%p)"
           "\n bad input\n", dv, fp) ;
   exit(-1) ;
}
size = DV_size(dv) ;
rc = fwrite((void *) &size, sizeof(int), 1, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in DV_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", dv, fp, rc, 1) ;
   return(0) ;
}
rc = fwrite((void *) DV_entries(dv), sizeof(double), size, fp) ;
if ( rc != size ) {
   fprintf(stderr, "\n error in DV_writeToBinaryFile(%p,%p)"
           "\n %d of %d items written\n",
           dv, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an DV object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   -------------------------------------------------
*/
int
DV_writeForHumanEye ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc ;

if ( dv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in DV_writeForHumanEye(%p,%p)"
           "\n bad input\n", dv, fp) ;
   exit(-1) ;
}
if ( (rc = DV_writeStats(dv, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in DV_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from DV_writeStats(%p,%p)\n",
           dv, fp, rc, dv, fp) ;
   return(0) ;
}
DVfprintf(fp, DV_size(dv), DV_entries(dv)) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the DV object

   return value -- 1 if success, 0 otherwise

   created -- 96jun23, cca
   ---------------------------------------------------------
*/
int
DV_writeStats ( 
   DV    *dv, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DV_writeStats(%p,%p)"
           "\n bad input\n", dv, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n DV : double vector object : ") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, 
             " size = %d, maxsize = %d, owned = %d", 
             dv->size, dv->maxsize, dv->owned) ;
if ( rc < 0 ) { goto IO_error ; }
return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in DV_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", dv, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write the DV object for a matlab file

   return value -- 1 if success, 0 otherwise

   created -- 98feb07, cca
   ---------------------------------------------------
*/
int
DV_writeForMatlab ( 
   DV     *dv, 
   char   *name,
   FILE   *fp 
) {
double   *entries ;
int      ii, rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DV_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", dv, name, fp) ;
   exit(-1) ;
}
DV_sizeAndEntries(dv, &size, &entries) ;
fprintf(fp, "\n %s = zeros(%d,%d) ;", name, size, size) ;
for ( ii = 0 ; ii < size ; ii++ ) {
   fprintf(fp, "\n %s(%d) = %24.16e ;", name, ii+1, entries[ii]) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
