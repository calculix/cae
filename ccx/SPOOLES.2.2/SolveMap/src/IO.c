/*  IO.c  */

#include "../SolveMap.h"

static const char *suffixb = ".solvemapb" ;
static const char *suffixf = ".solvemapf" ;

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to read in an SolveMap object from a file

   input --

      fn -- filename, must be *.solvemapb or *.solvemapf

   return value -- 1 if success, 0 if failure

   created -- 98apr09, cca
   -----------------------------------------------------
*/
int
SolveMap_readFromFile ( 
   SolveMap   *solvemap, 
   char       *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in SolveMap_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", solvemap, fn, __FILE__, __LINE__) ;
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
         fprintf(stderr, "\n error in SolveMap_readFromFile(%p,%s)"
                 "\n unable to open file %s", solvemap, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SolveMap_readFromBinaryFile(solvemap, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in SolveMap_readFromFile(%p,%s)"
                 "\n unable to open file %s", solvemap, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SolveMap_readFromFormattedFile(solvemap, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in SolveMap_readFromFile(%p,%s)"
              "\n bad SolveMap file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              solvemap, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in SolveMap_readFromFile(%p,%s)"
       "\n bad SolveMap file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       solvemap, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to read an SolveMap object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 98apr09, cca
   -----------------------------------------------------------
*/
int
SolveMap_readFromFormattedFile ( 
   SolveMap   *solvemap, 
   FILE       *fp 
) {
int   nblockLower, nblockUpper, nfront, nproc, rc, symmetryflag ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in SolveMap_readFromFormattedFile(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
SolveMap_clearData(solvemap) ;
/*
   -----------------------------------------------------
   read in the five scalar parameters
   symmetryflag, nfront, nproc, nblockUpper, nblockLower
   -----------------------------------------------------
*/
if ( (rc = IVfscanf(fp, 5, itemp)) != 5 ) {
   fprintf(stderr, "\n error in SolveMap_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", solvemap, fp, rc, 5) ;
   return(0) ;
}
symmetryflag = itemp[0] ;
nfront       = itemp[1] ;
nproc        = itemp[2] ;
nblockUpper  = itemp[3] ;
nblockLower  = itemp[4] ;
/*
   ---------------------
   initialize the object
   ---------------------
*/
SolveMap_init(solvemap, symmetryflag, nfront, nproc,
              nblockUpper, nblockLower) ;
/*
   ---------------------------
   read in the owners[] vector
   ---------------------------
*/
if ( (rc = IVfscanf(fp, nfront, solvemap->owners)) != nfront ) {
   fprintf(stderr, "\n error in SolveMap_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", solvemap, fp, rc, nfront) ;
   return(0) ;
}
if ( nblockUpper > 0 ) {
/*
   --------------------------------
   read in the rowidsUpper[] vector
   --------------------------------
*/
   if ( (rc = IVfscanf(fp, nblockUpper, solvemap->rowidsUpper)) 
       != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
/*
   --------------------------------
   read in the colidsUpper[] vector
   --------------------------------
*/
   if ( (rc = IVfscanf(fp, nblockUpper, solvemap->colidsUpper)) 
       != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
/*
   -----------------------------
   read in the mapUpper[] vector
   -----------------------------
*/
   if ( (rc = IVfscanf(fp, nblockUpper, solvemap->mapUpper)) 
       != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
}
if ( symmetryflag == 2 && nblockLower > 0 ) {
/*
   --------------------------------
   read in the rowidsLower[] vector
   --------------------------------
*/
   if ( (rc = IVfscanf(fp, nblockLower, solvemap->rowidsLower)) 
       != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
/*
   --------------------------------
   read in the colidsLower[] vector
   --------------------------------
*/
   if ( (rc = IVfscanf(fp, nblockLower, solvemap->colidsLower)) 
       != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
/*
   -----------------------------
   read in the mapLower[] vector
   -----------------------------
*/
   if ( (rc = IVfscanf(fp, nblockLower, solvemap->mapLower)) 
       != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to read an SolveMap object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 98apr09, cca
   ---------------------------------------------------------
*/
int
SolveMap_readFromBinaryFile ( 
   SolveMap    *solvemap, 
   FILE   *fp 
) {
int   nblockLower, nblockUpper, nfront, nproc, rc, symmetryflag ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_readFromBinaryFile(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
SolveMap_clearData(solvemap) ;
/*
   -----------------------------------------------------
   read in the five scalar parameters
   symmetryflag, nfront, nproc, nblockUpper, nblockLower
   -----------------------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 5, fp)) != 5 ) {
   fprintf(stderr, "\n error in SolveMap_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", solvemap, fp, rc, 5) ;
   return(0) ;
}
symmetryflag = itemp[0] ;
nfront       = itemp[1] ;
nproc        = itemp[2] ;
nblockUpper  = itemp[3] ;
nblockLower  = itemp[4] ;
/*
   ---------------------
   initialize the object
   ---------------------
*/
SolveMap_init(solvemap, symmetryflag, nfront, nproc, 
              nblockUpper, nblockLower) ;
/*
   ---------------------------
   read in the owners[] vector
   ---------------------------
*/
fread((void *) solvemap->owners, sizeof(int), nfront, fp) ;
if ( rc != nfront ) {
   fprintf(stderr, "\n error in SolveMap_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", solvemap, fp, rc, nfront) ;
   return(0) ;
}
if ( nblockUpper > 0 ) {
/*
   --------------------------------
   read in the rowidsUpper[] vector
   --------------------------------
*/
   fread((void *) solvemap->rowidsUpper, sizeof(int), nblockUpper, fp) ;
   if ( rc != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
/*
   --------------------------------
   read in the colidsUpper[] vector
   --------------------------------
*/
   fread((void *) solvemap->colidsUpper, sizeof(int), nblockUpper, fp) ;
   if ( rc != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
/*
   -----------------------------
   read in the mapUpper[] vector
   -----------------------------
*/
   fread((void *) solvemap->mapUpper, sizeof(int), nblockUpper, fp) ;
   if ( rc != nblockUpper ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockUpper) ;
      return(0) ;
   }
}
if ( symmetryflag == 2 && nblockLower > 0 ) {
/*
   --------------------------------
   read in the rowidsLower[] vector
   --------------------------------
*/
   fread((void *) solvemap->rowidsLower, sizeof(int), nblockLower, fp) ;
   if ( rc != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
/*
   --------------------------------
   read in the colidsLower[] vector
   --------------------------------
*/
   fread((void *) solvemap->colidsLower, sizeof(int), nblockLower, fp) ;
   if ( rc != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
/*
   -----------------------------
   read in the mapLower[] vector
   -----------------------------
*/
   fread((void *) solvemap->mapLower, sizeof(int), nblockLower, fp) ;
   if ( rc != nblockLower ) {
      fprintf(stderr, 
              "\n error in SolveMap_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", 
              solvemap, fp, rc, nblockLower) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to write an SolveMap object to a file

   input --

      fn -- filename
        *.solvemapb -- binary
        *.solvemapf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 98apr09, cca
   ------------------------------------------------
*/
int
SolveMap_writeToFile ( 
   SolveMap    *solvemap, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_writeToFile(%p,%s)"
    "\n bad input\n", solvemap, fn) ; 
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
         fprintf(stderr, "\n error in SolveMap_writeToFile(%p,%s)"
                 "\n unable to open file %s", solvemap, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SolveMap_writeToBinaryFile(solvemap, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in SolveMap_writeToFile(%p,%s)"
                 "\n unable to open file %s", solvemap, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SolveMap_writeToFormattedFile(solvemap, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in SolveMap_writeToFile(%p,%s)"
                 "\n unable to open file %s", solvemap, fn, fn) ;
         rc = 0 ;
      } else {
         rc = SolveMap_writeForHumanEye(solvemap, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in SolveMap_writeToFile(%p,%s)"
              "\n unable to open file %s", solvemap, fn, fn) ;
      rc = 0 ;
   } else {
      rc = SolveMap_writeForHumanEye(solvemap, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an SolveMap object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 98apr09, cca
   -----------------------------------------------------
*/
int
SolveMap_writeToFormattedFile ( 
   SolveMap    *solvemap, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   exit(-1) ;
}
/*
   ------------------------------------
   write out the five scalar parameters
   ------------------------------------
*/
rc = fprintf(fp, "\n %d %d %d %d %d", 
             solvemap->symmetryflag, 
             solvemap->nfront, 
             solvemap->nproc, 
             solvemap->nblockUpper, 
             solvemap->nblockLower) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", 
           solvemap, fp, rc) ;
   return(0) ;
}
if ( solvemap->nfront > 0 ) {
   IVfp80(fp, solvemap->nfront, solvemap->owners, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from owners[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
}
if ( solvemap->nblockUpper > 0 ) {
   IVfp80(fp, solvemap->nblockUpper, solvemap->rowidsUpper, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from rowidsUpper[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
   IVfp80(fp, solvemap->nblockUpper, solvemap->colidsUpper, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from colidsUpper[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
   IVfp80(fp, solvemap->nblockUpper, solvemap->mapUpper, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from mapUpper[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
}
if ( solvemap->symmetryflag == 2 && solvemap->nblockLower > 0 ) {
   IVfp80(fp, solvemap->nblockLower, solvemap->rowidsLower, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from rowidsLower[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
   IVfp80(fp, solvemap->nblockLower, solvemap->colidsLower, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from colidsLower[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
   IVfp80(fp, solvemap->nblockLower, solvemap->mapLower, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in SolveMap_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from mapLower[] IVfp80\n", 
              solvemap, fp, ierr) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an SolveMap object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 98apr09, cca
   --------------------------------------------------
*/
int
SolveMap_writeToBinaryFile ( 
   SolveMap    *solvemap, 
   FILE   *fp 
) {
int   rc ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_writeToBinaryFile(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   exit(-1) ;
}
itemp[0] = solvemap->symmetryflag ;
itemp[1] = solvemap->nfront       ;
itemp[2] = solvemap->nproc        ;
itemp[2] = solvemap->nblockUpper  ;
itemp[2] = solvemap->nblockLower  ;
rc = fwrite((void *) itemp, sizeof(int), 5, fp) ;
if ( rc != 5 ) {
   fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", solvemap, fp, rc, 5) ;
   return(0) ;
}
rc = fwrite((void *) solvemap->owners, sizeof(int), 
             solvemap->nfront, fp) ;
if ( rc != solvemap->nfront ) {
   fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
           "\n owners: %d of %d items written\n",
           solvemap, fp, rc, solvemap->nfront) ;
   return(0) ;
}
if ( solvemap->nblockUpper > 0 ) {
   rc = fwrite((void *) solvemap->rowidsUpper, sizeof(int), 
                solvemap->nblockUpper, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockUpper) ;
      return(0) ;
   }
   rc = fwrite((void *) solvemap->colidsUpper, sizeof(int), 
                solvemap->nblockUpper, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockUpper) ;
      return(0) ;
   }
   rc = fwrite((void *) solvemap->mapUpper, sizeof(int), 
                solvemap->nblockUpper, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockUpper) ;
      return(0) ;
   }
}
if ( solvemap->symmetryflag == 2 && solvemap->nblockLower > 0 ) {
   rc = fwrite((void *) solvemap->rowidsLower, sizeof(int), 
                solvemap->nblockLower, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockLower) ;
      return(0) ;
   }
   rc = fwrite((void *) solvemap->colidsLower, sizeof(int), 
                solvemap->nblockLower, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockLower) ;
      return(0) ;
   }
   rc = fwrite((void *) solvemap->mapLower, sizeof(int), 
                solvemap->nblockLower, fp) ;
   if ( rc != solvemap->nfront ) {
      fprintf(stderr, "\n error in SolveMap_writeToBinaryFile(%p,%p)"
              "\n owners: %d of %d items written\n",
              solvemap, fp, rc, solvemap->nblockLower) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write an SolveMap object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 98apr09, cca
   ------------------------------------------------------
*/
int
SolveMap_writeForHumanEye ( 
   SolveMap   *solvemap, 
   FILE       *fp 
) {
int   ierr, kk, rc ;

if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in SolveMap_writeForHumanEye(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   exit(-1) ;
}
if ( (rc = SolveMap_writeStats(solvemap, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in SolveMap_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from SolveMap_writeStats(%p,%p)\n",
           solvemap, fp, rc, solvemap, fp) ;
   return(0) ;
}
fprintf(fp, "\n\n front owners map") ;
IVfp80(fp, solvemap->nfront, solvemap->owners, 80, &ierr) ;
if ( solvemap->nblockUpper > 0 ) {
   fprintf(fp, "\n\n upper block map") ;
   for ( kk = 0 ; kk < solvemap->nblockUpper ; kk++ ) {
      fprintf(fp, "\n block(%d,%d) owned by process %d",
              solvemap->rowidsUpper[kk], solvemap->colidsUpper[kk],
              solvemap->mapUpper[kk]) ;
   }
}
if ( solvemap->symmetryflag == 2 && solvemap->nblockLower > 0 ) {
   fprintf(fp, "\n\n lower block map") ;
   for ( kk = 0 ; kk < solvemap->nblockLower ; kk++ ) {
      fprintf(fp, "\n block(%d,%d) owned by process %d",
              solvemap->rowidsLower[kk], solvemap->colidsLower[kk],
              solvemap->mapLower[kk]) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- to write out the statistics for the SolveMap object

   return value -- 1 if success, 0 otherwise

   created -- 98apr09, cca
   --------------------------------------------------------------
*/
int
SolveMap_writeStats ( 
   SolveMap   *solvemap, 
   FILE       *fp 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( solvemap == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in SolveMap_writeStats(%p,%p)"
           "\n bad input\n", solvemap, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n SolveMap : submatrix solve object :") ;
if ( solvemap->symmetryflag < 2 ) {
   fprintf(fp, "\n matrix is symmetric") ;
} else {
   fprintf(fp, "\n matrix is nonsymmetric") ;
}
fprintf(fp, 
         "\n %d fronts, %d processes, %d upper blocks, %d lower blocks",
         solvemap->nfront, solvemap->nproc, 
         solvemap->nblockUpper, solvemap->nblockLower) ;
return(1) ; }

/*--------------------------------------------------------------------*/
