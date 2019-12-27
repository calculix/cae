/*  IO.c  */

#include "../Perm.h"

static const char *suffixb = ".permb" ;
static const char *suffixf = ".permf" ;

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to read in an Perm object from a file

   input --

      fn -- filename, must be *.permb or *.permf

   return value -- 1 if success, 0 if failure

   created -- 96jan05, cca
   -----------------------------------------------
*/
int
Perm_readFromFile ( 
   Perm   *perm, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in Perm_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", perm, fn, __FILE__, __LINE__) ;
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
         fprintf(stderr, "\n error in Perm_readFromFile(%p,%s)"
                 "\n unable to open file %s", perm, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Perm_readFromBinaryFile(perm, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in Perm_readFromFile(%p,%s)"
                 "\n unable to open file %s", perm, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Perm_readFromFormattedFile(perm, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in Perm_readFromFile(%p,%s)"
              "\n bad Perm file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              perm, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in Perm_readFromFile(%p,%s)"
       "\n bad Perm file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       perm, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read an Perm object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96jan05, cca
   ------------------------------------------------------
*/
int
Perm_readFromFormattedFile ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   i, isPresent, j, rc, size ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Perm_readFromFormattedFile(%p,%p)"
           "\n bad input\n", perm, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Perm_clearData(perm) ;
/*
   -----------------------------------------------------
   read in the two scalar parameters: isPresent and size
   -----------------------------------------------------
*/
if ( (rc = IVfscanf(fp, 2, itemp)) != 2 ) {
   fprintf(stderr, "\n error in Perm_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", perm, fp, rc, 2) ;
   return(0) ;
}
isPresent = itemp[0] ;
size      = itemp[1] ;
if ( isPresent < 1 || isPresent > 3 || size <= 0 ) {
   fprintf(stderr, "\n error in Perm_readFromFormattedFile(%p,%p)"
           "\n isPresent = %d, size = %d", perm, fp, isPresent, size) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Perm_initWithTypeAndSize(perm, isPresent, size) ;
if ( isPresent == 2 || isPresent == 3 ) {
/*
   -----------------------------------------
   read in the old-to-new permutation vector
   -----------------------------------------
*/
   int   *oldToNew = perm->oldToNew ;
   if ( (rc = IVfscanf(fp, size, oldToNew)) != size ) {
      fprintf(stderr, "\n error in Perm_readFromFormattedFile(%p,%p)"
             "\n %d items of %d read\n", 
             perm, fp, rc, size) ;
      exit(-1) ;
   }
   for ( i = 0 ; i < size ; i++ ) {
      if ( oldToNew[i] == size ) {
         for ( j = 0 ; j < size ; j++ ) {
            oldToNew[j]-- ;
         }
         break ;
      }
   }
}
if ( isPresent == 1 || isPresent == 3 ) {
/*
   -----------------------------------------
   read in the new-to-old permutation vector
   -----------------------------------------
*/
   int   *newToOld = perm->newToOld ;
   if ( (rc = IVfscanf(fp, size, newToOld)) != size ) {
      fprintf(stderr, "\n error in Perm_readFromFormattedFile(%p,%p)"
             "\n %d items of %d read\n", 
             perm, fp, rc, size) ;
      exit(-1) ;
   }
   for ( i = 0 ; i < size ; i++ ) {
      if ( newToOld[i] == size ) {
         for ( j = 0 ; j < size ; j++ ) {
            newToOld[j]-- ;
         }
         break ;
      }
   }
}
/*
   ----------------------------
   check the permutation object
   ----------------------------
*/
if ( Perm_checkPerm(perm) != 1 ) {
   fprintf(stderr, "\n fatal error in Perm_readFromFormattedFile(%p,%p)"
           "\n permutation is not valid\n", perm, fp) ;
   exit(-1) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read an Perm object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96jan05, cca
   ---------------------------------------------------
*/
int
Perm_readFromBinaryFile ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   i, isPresent, j, rc, size ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_readFromBinaryFile(%p,%p)"
           "\n bad input\n", perm, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Perm_clearData(perm) ;
/*
   ------------------------------------------------
   read in the two scalar parameters: type and size
   ------------------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 2, fp)) != 2 ) {
   fprintf(stderr, "\n error in Perm_readFromBinaryFile(%p,%p)"
           "\n itemp(2) : %d items of %d read\n", perm, fp, rc, 2) ;
   return(0) ;
}
isPresent = itemp[0] ;
size      = itemp[1] ;
if ( isPresent < 1 || isPresent > 3 || size <= 0 ) {
   fprintf(stderr, "\n error in Perm_readFromBinaryFile(%p,%p)"
           "\n isPresent = %d, size = %d", perm, fp, isPresent, size) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Perm_initWithTypeAndSize(perm, isPresent, size) ;
if ( isPresent == 2 || isPresent == 3 ) {
/*
   -----------------------------------------
   read in the old-to-new permutation vector
   -----------------------------------------
*/
   int   *oldToNew = perm->oldToNew ;
   if ( (rc = fread(oldToNew, sizeof(int), size, fp)) != size ) {
      fprintf(stderr, "\n error in Perm_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", perm, fp, rc, size) ;
      exit(-1) ;
   }
   for ( i = 0 ; i < size ; i++ ) {
      if ( oldToNew[i] == size ) {
         for ( j = 0 ; j < size ; j++ ) {
            oldToNew[j]-- ;
         }
         break ;
      }
   }
}
if ( isPresent == 1 || isPresent == 3 ) {
/*
   -----------------------------------------
   read in the new-to-old permutation vector
   -----------------------------------------
*/
   int   *newToOld = perm->newToOld ;
   if ( (rc = fread(newToOld, sizeof(int), size, fp)) != size ) {
      fprintf(stderr, "\n error in Perm_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", perm, fp, rc, size) ;
      exit(-1) ;
   }
   for ( i = 0 ; i < size ; i++ ) {
      if ( newToOld[i] == size ) {
         for ( j = 0 ; j < size ; j++ ) {
            newToOld[j]-- ;
         }
         break ;
      }
   }
}
if ( Perm_checkPerm(perm) != 1 ) {
   fprintf(stderr, "\n fatal error in Perm_readFromFormattedFile(%p,%p)"
           "\n permutation is not valid\n", perm, fp) ;
   exit(-1) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write an Perm object to a file

   input --

      fn -- filename
        *.permb -- binary
        *.permf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -------------------------------------------
*/
int
Perm_writeToFile ( 
   Perm   *perm, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_writeToFile(%p,%s)"
    "\n bad input\n", perm, fn) ; 
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
         fprintf(stderr, "\n error in Perm_writeToFile(%p,%s)"
                 "\n unable to open file %s", perm, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Perm_writeToBinaryFile(perm, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in Perm_writeToFile(%p,%s)"
                 "\n unable to open file %s", perm, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Perm_writeToFormattedFile(perm, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in Perm_writeToFile(%p,%s)"
                 "\n unable to open file %s", perm, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Perm_writeForHumanEye(perm, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in Perm_writeToFile(%p,%s)"
              "\n unable to open file %s", perm, fn, fn) ;
      rc = 0 ;
   } else {
      rc = Perm_writeForHumanEye(perm, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an Perm object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -----------------------------------------------------
*/
int
Perm_writeToFormattedFile ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fp == NULL || perm->size <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_writeToFormattedFile(%p,%p)"
           "\n bad input\n", perm, fp) ;
   exit(-1) ;
}
/*
   -----------------------------------
   write out the two scalar parameters
   -----------------------------------
*/
rc = fprintf(fp, "\n %d %d", perm->isPresent, perm->size) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in Perm_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", perm, fp, rc) ;
   return(0) ;
}
if ( perm->isPresent == 2 || perm->isPresent == 3 ) {
/*
   -------------------------------------------
   write out the old-to-new permutation vector
   -------------------------------------------
*/
   IVfp80(fp, perm->size, perm->oldToNew, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Perm_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from oldToNew[] IVfp80\n", 
              perm, fp, ierr) ;
      return(0) ;
   }
}
if ( perm->isPresent == 1 || perm->isPresent == 3 ) {
/*
   -------------------------------------------
   write out the new-to-old permutation vector
   -------------------------------------------
*/
   IVfp80(fp, perm->size, perm->newToOld, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Perm_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from newToOld[] IVfp80\n", 
              perm, fp, ierr) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an Perm object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   --------------------------------------------------
*/
int
Perm_writeToBinaryFile ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   size, rc ;
int   itemp[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fp == NULL || (size = perm->size) <= 0 ) {
   fprintf(stderr, "\n fatal error in Perm_writeToBinaryFile(%p,%p)"
           "\n bad input\n", perm, fp) ;
   exit(-1) ;
}
/*
   -----------------------------------
   write out the two scalar parameters
   -----------------------------------
*/
itemp[0] = perm->isPresent ;
itemp[1] = size ;
rc = fwrite((void *) itemp, sizeof(int), 2, fp) ;
if ( rc != 2 ) {
   fprintf(stderr, "\n error in Perm_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", perm, fp, rc, 2) ;
   return(0) ;
}
if ( perm->isPresent == 2 || perm->isPresent == 3 ) {
/*
   -------------------------------------------
   write out the old-to-new permutation vector
   -------------------------------------------
*/
   rc = fwrite((void *) perm->oldToNew, sizeof(int), size, fp) ;
   if ( rc != size ) {
      fprintf(stderr, "\n error in Perm_writeToBinaryFile(%p,%p)"
              "\n perm->oldToNew, %d of %d items written\n",
              perm, fp, rc, size) ;
      return(0) ;
   }
}
if ( perm->isPresent == 1 || perm->isPresent == 3 ) {
/*
   -------------------------------------------
   write out the new-to-old permutation vector
   -------------------------------------------
*/
   rc = fwrite((void *) perm->newToOld, sizeof(int), size, fp) ;
   if ( rc != size ) {
      fprintf(stderr, "\n error in Perm_writeToBinaryFile(%p,%p)"
              "\n perm->newToOld, %d of %d items written\n",
              perm, fp, rc, size) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an Perm object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   -------------------------------------------------
*/
int
Perm_writeForHumanEye ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   ierr, rc ;

if ( perm == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Perm_writeForHumanEye(%p,%p)"
           "\n bad input\n", perm, fp) ;
   exit(-1) ;
}
if ( (rc = Perm_writeStats(perm, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in Perm_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from Perm_writeStats(%p,%p)\n",
           perm, fp, rc, perm, fp) ;
   return(0) ;
}
if ( perm->isPresent == 2 || perm->isPresent == 3 ) {
   fprintf(fp, "\n\n old-to-new permutation") ;
   IVfp80(fp, perm->size, perm->oldToNew, 80, &ierr) ;
}
if ( perm->isPresent == 1 || perm->isPresent == 3 ) {
   fprintf(fp, "\n\n new-to-old permutation") ;
   IVfp80(fp, perm->size, perm->newToOld, 80, &ierr) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the Perm object

   return value -- 1 if success, 0 otherwise

   created -- 96jan05, cca
   ---------------------------------------------------------
*/
int
Perm_writeStats ( 
   Perm   *perm, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( perm == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Perm_writeStats(%p,%p)"
           "\n bad input\n", perm, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n Perm : permutation object :") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, "\n isPresent %d, size %d", 
             perm->isPresent, perm->size) ;
if ( rc < 0 ) { goto IO_error ; }
return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in Perm_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", perm, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
