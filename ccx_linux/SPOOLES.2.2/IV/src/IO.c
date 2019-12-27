/*  IO.c  */

#include "../IV.h"

static const char *suffixb = ".ivb" ;
static const char *suffixf = ".ivf" ;

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to read in an IV object from a file

   input --

      fn -- filename, must be *.ivb or *.ivf

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
IV_readFromFile ( 
   IV    *iv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in IV_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", iv, fn, __FILE__, __LINE__) ;
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
         fprintf(stderr, "\n error in IV_readFromFile(%p,%s)"
                 "\n unable to open file %s", iv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IV_readFromBinaryFile(iv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in IV_readFromFile(%p,%s)"
                 "\n unable to open file %s", iv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IV_readFromFormattedFile(iv, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in IV_readFromFile(%p,%s)"
              "\n bad IV file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              iv, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in IV_readFromFile(%p,%s)"
       "\n bad IV file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       iv, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to read an IV object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   -----------------------------------------------------
*/
int
IV_readFromFormattedFile ( 
   IV    *iv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in IV_readFromFormattedFile(%p,%p)"
           "\n bad input\n", iv, fp) ;
   return(0) ;
}
IV_clearData(iv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fscanf(fp, "%d", &size)) != 1 ) {
   fprintf(stderr, "\n error in IV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", iv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
IV_init(iv, size, NULL) ;
iv->size = size ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = IVfscanf(fp, size, iv->vec)) != size ) {
   fprintf(stderr, "\n error in IV_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", iv, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read an IV object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95oct06, cca
   ---------------------------------------------------
*/
int
IV_readFromBinaryFile ( 
   IV    *iv, 
   FILE   *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in IV_readFromBinaryFile(%p,%p)"
           "\n bad input\n", iv, fp) ;
   return(0) ;
}
IV_clearData(iv) ;
/*
   ------------------------------
   read in the size of the vector
   ------------------------------
*/
if ( (rc = fread((void *) &size, sizeof(int), 1, fp)) != 1 ) {
   fprintf(stderr, "\n error in IV_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", iv, fp, rc, 1) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
IV_init(iv, size, NULL) ;
iv->size = size ;
/*
   ------------------------
   read in the vec[] vector
   ------------------------
*/
if ( (rc = fread((void *) iv->vec, sizeof(int), size, fp)) != size ) {
   fprintf(stderr, "\n error in IV_readFromBinaryFile(%p,%p)"
           "\n sizes(%d) : %d items of %d read\n", 
           iv, fp, size, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write an IV object to a file

   input --

      fn -- filename
        *.ivb -- binary
        *.ivf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------
*/
int
IV_writeToFile ( 
   IV    *iv, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in IV_writeToFile(%p,%s)"
    "\n bad input\n", iv, fn) ; 
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
         fprintf(stderr, "\n error in IV_writeToFile(%p,%s)"
                 "\n unable to open file %s", iv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IV_writeToBinaryFile(iv, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in IV_writeToFile(%p,%s)"
                 "\n unable to open file %s", iv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IV_writeToFormattedFile(iv, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in IV_writeToFile(%p,%s)"
                 "\n unable to open file %s", iv, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IV_writeForHumanEye(iv, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in IV_writeToFile(%p,%s)"
              "\n unable to open file %s", iv, fn, fn) ;
      rc = 0 ;
   } else {
      rc = IV_writeForHumanEye(iv, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an IV object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -----------------------------------------------------
*/
int
IV_writeToFormattedFile ( 
   IV    *iv, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL || iv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in IV_writeToFormattedFile(%p,%p)"
           "\n bad input\n", iv, fp) ;
   fprintf(stderr, "\n iv->size = %d", iv->size) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the size of the vector
   -------------------------------------
*/
rc = fprintf(fp, "\n %d", iv->size) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in IV_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", iv, fp, rc) ;
   return(0) ;
}
if ( iv->size > 0 ) {
   IVfp80(fp, iv->size, iv->vec, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in IV_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from sizes[] IVfp80\n", 
              iv, fp, ierr) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an IV object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   --------------------------------------------------
*/
int
IV_writeToBinaryFile ( 
   IV    *iv, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL || iv->size <= 0 ) {
   fprintf(stderr, "\n fatal error in IV_writeToBinaryFile(%p,%p)"
           "\n bad input\n", iv, fp) ;
   exit(-1) ;
}
rc = fwrite((void *) &iv->size, sizeof(int), 1, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in IV_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", iv, fp, rc, 1) ;
   return(0) ;
}
rc = fwrite((void *) iv->vec, sizeof(int), iv->size, fp) ;
if ( rc != iv->size ) {
   fprintf(stderr, "\n error in IV_writeToBinaryFile(%p,%p)"
           "\n iv->sizes, %d of %d items written\n",
           iv, fp, rc, iv->size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an IV object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------------
*/
int
IV_writeForHumanEye ( 
   IV    *iv, 
   FILE   *fp 
) {
int   ierr, rc ;

if ( iv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in IV_writeForHumanEye(%p,%p)"
           "\n bad input\n", iv, fp) ;
   exit(-1) ;
}
if ( (rc = IV_writeStats(iv, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in IV_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from IV_writeStats(%p,%p)\n",
           iv, fp, rc, iv, fp) ;
   return(0) ;
}
IVfp80(fp, iv->size, iv->vec, 80, &ierr) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the IV object

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ---------------------------------------------------------
*/
int
IV_writeStats ( 
   IV    *iv, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in IV_writeStats(%p,%p)"
           "\n bad input\n", iv, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n IV : integer vector object : ") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, " size = %d, maxsize = %d, owned = %d", 
             iv->size, iv->maxsize, iv->owned) ;
if ( rc < 0 ) { goto IO_error ; }
return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in IV_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", iv, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
  purpose -- to write out an integer vector with eighty column lines
 
   input --
 
      fp     -- file pointer, must be formatted and write access
      column -- present column
      pierr  -- pointer to int to hold return value, 
                should be 1 if any print was successful,
                if fprintf() failed, then ierr = -1
  
   return value -- present column
 
   created -- 96jun22, cca
   -------------------------------------------------------------------
*/
int
IV_fp80 (
   IV     *iv,
   FILE   *fp,
   int    column,
   int    *pierr
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL || pierr == NULL ) {
   fprintf(stderr, "\n fatal error in IV_fp80(%p,%p,%p)"
           "\n bad input\n", iv, fp, pierr) ;
   exit(-1) ;
}
if ( iv->size > 0 && iv->vec != NULL ) {
   column = IVfp80(fp, iv->size, iv->vec, column, pierr) ;
}

return(column) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write the IV object for a matlab file
 
   return value -- 1 if success, 0 otherwise
 
   created -- 98oct22, cca
   ---------------------------------------------------
*/
int
IV_writeForMatlab ( 
   IV     *iv, 
   char   *name,
   FILE   *fp 
) {
int   ii, rc, size ;
int   *entries ;
/*
   ---------------
   check the input
   ---------------
*/
if ( iv == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in IV_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", iv, name, fp) ;
   exit(-1) ;
}
IV_sizeAndEntries(iv, &size, &entries) ;
fprintf(fp, "\n %s = zeros(%d,%d) ;", name, size, size) ;
for ( ii = 0 ; ii < size ; ii++ ) {
   fprintf(fp, "\n %s(%d) = %d ;", name, ii+1, entries[ii]+1) ;
}
return(1) ; }
 
/*--------------------------------------------------------------------*/
