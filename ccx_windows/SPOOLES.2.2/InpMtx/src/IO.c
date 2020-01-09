/*  IO.c  */

#include "../../InpMtx.h"

static const char *suffixb = ".inpmtxb" ;
static const char *suffixf = ".inpmtxf" ;

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read in a InpMtx object from a file

   input --

      fn -- filename, must be *.inpmtxb or *.inpmtxf

   return value -- 1 if success, 0 if failure

   created -- 98jan28, cca
   --------------------------------------------------
*/
int
InpMtx_readFromFile ( 
   InpMtx   *inpmtx, 
   char      *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
           "\n bad input\n", inpmtx, fn) ;
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
         fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", inpmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = InpMtx_readFromBinaryFile(inpmtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", inpmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = InpMtx_readFromFormattedFile(inpmtx, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
              "\n bad InpMtx file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              inpmtx, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
       "\n bad InpMtx file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       inpmtx, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read a InpMtx object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 98jan28, cca
   --------------------------------------------------------
*/
int
InpMtx_readFromFormattedFile ( 
   InpMtx   *inpmtx, 
   FILE      *fp 
) {
int   rc ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in InpMtx_readFromFormattedFile(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
InpMtx_clearData(inpmtx) ;
/*
   --------------------------------------------------------
   read in the five scalar parameters
   coordinate type, storage mode, input mode, nent, nvector
   --------------------------------------------------------
*/
if ( (rc = IVfscanf(fp, 5, itemp)) != 5 ) {
   fprintf(stderr, "\n error in InpMtx_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", inpmtx, fp, rc, 5) ;
   return(0) ;
}
inpmtx->coordType   = itemp[0] ;
inpmtx->storageMode = itemp[1] ;
inpmtx->inputMode   = itemp[2] ;
inpmtx->nent        = itemp[3] ;
inpmtx->nvector     = itemp[4] ;
if ( inpmtx->nent > 0 ) {
   IV_readFromFormattedFile(&inpmtx->ivec1IV, fp) ;
   IV_readFromFormattedFile(&inpmtx->ivec2IV, fp) ;
   if (  INPMTX_IS_REAL_ENTRIES(inpmtx) 
      || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      DV_readFromFormattedFile(&inpmtx->dvecDV, fp) ;
   }
}
if ( inpmtx->nvector > 0 ) {
   IV_readFromFormattedFile(&inpmtx->vecidsIV,  fp) ;
   IV_readFromFormattedFile(&inpmtx->sizesIV,   fp) ;
   IV_readFromFormattedFile(&inpmtx->offsetsIV, fp) ;
}
inpmtx->maxnent = inpmtx->nent ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read a InpMtx object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 98jan28, cca
   ------------------------------------------------------
*/
int
InpMtx_readFromBinaryFile ( 
   InpMtx   *inpmtx, 
   FILE    *fp 
) {
int   rc ;
int   itemp[5] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_readFromBinaryFile(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
InpMtx_clearData(inpmtx) ;
/*
   ---------------------------------------------
   read in the five scalar parameters
   coordType storageMode inputMode nent nvector
   ---------------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 5, fp)) != 5 ) {
   fprintf(stderr, "\n error in InpMtx_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", inpmtx, fp, rc, 5) ;
   return(0) ;
}
inpmtx->coordType   = itemp[0] ;
inpmtx->storageMode = itemp[1] ;
inpmtx->inputMode   = itemp[2] ;
inpmtx->nent        = itemp[3] ;
inpmtx->nvector     = itemp[4] ;
if ( inpmtx->nent > 0 ) {
   IV_readFromBinaryFile(&inpmtx->ivec1IV, fp) ;
   IV_readFromBinaryFile(&inpmtx->ivec2IV, fp) ;
   if (  INPMTX_IS_REAL_ENTRIES(inpmtx) 
      || INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      DV_readFromBinaryFile(&inpmtx->dvecDV, fp) ;
   }
}
if ( inpmtx->nvector > 0 ) {
   IV_readFromBinaryFile(&inpmtx->vecidsIV,  fp) ;
   IV_readFromBinaryFile(&inpmtx->sizesIV,   fp) ;
   IV_readFromBinaryFile(&inpmtx->offsetsIV, fp) ;
}
inpmtx->maxnent = inpmtx->nent ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to write a InpMtx object to a file

   input --

      fn -- filename
        *.inpmtxb -- binary
        *.inpmtxf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   ----------------------------------------------
*/
int
InpMtx_writeToFile ( 
   InpMtx   *inpmtx, 
   char    *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_writeToFile(%p,%s)"
    "\n bad input\n", inpmtx, fn) ; 
   return(0) ;
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
         fprintf(stderr, "\n error in InpMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", inpmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = InpMtx_writeToBinaryFile(inpmtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in InpMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", inpmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = InpMtx_writeToFormattedFile(inpmtx, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in InpMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", inpmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = InpMtx_writeForHumanEye(inpmtx, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in InpMtx_writeToFile(%p,%s)"
              "\n unable to open file %s", inpmtx, fn, fn) ;
      rc = 0 ;
   } else {
      rc = InpMtx_writeForHumanEye(inpmtx, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write a InpMtx object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   ------------------------------------------------------
*/
int
InpMtx_writeToFormattedFile ( 
   InpMtx   *inpmtx, 
   FILE    *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   return(0) ;
}
/*
   ------------------------------------
   write out the five scalar parameters
   ------------------------------------
*/
rc = fprintf(fp, "\n %d %d %d %d %d", 
             inpmtx->coordType, inpmtx->storageMode, 
             inpmtx->inputMode, inpmtx->nent, inpmtx->nvector) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", inpmtx, fp, rc) ;
   return(0) ;
}
/*
   ---------------------
   write out the vectors
   ---------------------
*/
if ( inpmtx->nent > 0 ) {
   IV_setSize(&inpmtx->ivec1IV, inpmtx->nent) ;
   rc = IV_writeToFormattedFile(&inpmtx->ivec1IV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from writing ivec1\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   IV_setSize(&inpmtx->ivec2IV, inpmtx->nent) ;
   rc = IV_writeToFormattedFile(&inpmtx->ivec2IV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from writing ivec2\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   if (  INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
      DV_setSize(&inpmtx->dvecDV, inpmtx->nent) ;
      rc = DV_writeToFormattedFile(&inpmtx->dvecDV,  fp) ;
      if ( rc < 0 ) {
         fprintf(stderr, 
                 "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
                 "\n rc = %d, return from writing dvec\n", 
                 inpmtx, fp, rc) ;
         return(0) ;
      }
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      DV_setSize(&inpmtx->dvecDV, 2*inpmtx->nent) ;
      rc = DV_writeToFormattedFile(&inpmtx->dvecDV,  fp) ;
      if ( rc < 0 ) {
         fprintf(stderr, 
                 "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
                 "\n rc = %d, return from writing dvec\n", 
                 inpmtx, fp, rc) ;
         return(0) ;
      }
   }
}
if ( inpmtx->nvector > 0 ) {
   rc = IV_writeToFormattedFile(&inpmtx->vecidsIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from writing vecids\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   rc = IV_writeToFormattedFile(&inpmtx->sizesIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from writing sizes\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   rc = IV_writeToFormattedFile(&inpmtx->offsetsIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from writing offsets\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write a InpMtx object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   ---------------------------------------------------
*/
int
InpMtx_writeToBinaryFile ( 
   InpMtx    *inpmtx, 
   FILE   *fp 
) {
int   rc ;
int   itemp[6] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   return(0) ;
}
/*
   ------------------------------------
   write out the five scalar parameters
   ------------------------------------
*/
itemp[0] = inpmtx->coordType   ;
itemp[1] = inpmtx->storageMode ;
itemp[2] = inpmtx->inputMode   ;
itemp[3] = inpmtx->nent        ;
itemp[4] = inpmtx->nvector     ;
rc = fwrite((void *) itemp, sizeof(int), 5, fp) ;
if ( rc != 5 ) {
   fprintf(stderr, "\n error in InpMtx_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", inpmtx, fp, rc, 5) ;
   return(0) ;
}
/*
   ---------------------
   write out the vectors
   ---------------------
*/
if ( inpmtx->nent > 0 ) {
   IV_setSize(&inpmtx->ivec1IV, inpmtx->nent) ;
   rc = IV_writeToBinaryFile(&inpmtx->ivec1IV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from writing ivec1\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   IV_setSize(&inpmtx->ivec2IV, inpmtx->nent) ;
   rc = IV_writeToBinaryFile(&inpmtx->ivec2IV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from writing ivec2\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
      DV_setSize(&inpmtx->dvecDV, inpmtx->nent) ;
      rc = DV_writeToBinaryFile(&inpmtx->dvecDV, fp) ;
      if ( rc < 0 ) {
         fprintf(stderr, 
                 "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
                 "\n rc = %d, return from writing dvec\n", 
                 inpmtx, fp, rc) ;
         return(0) ;
      }
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      DV_setSize(&inpmtx->dvecDV, 2*inpmtx->nent) ;
      rc = DV_writeToBinaryFile(&inpmtx->dvecDV, fp) ;
      if ( rc < 0 ) {
         fprintf(stderr, 
                 "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
                 "\n rc = %d, return from writing dvec\n", 
                 inpmtx, fp, rc) ;
         return(0) ;
      }
   }
}
if ( inpmtx->nvector > 0 ) {
   rc = IV_writeToBinaryFile(&inpmtx->vecidsIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from writing vecids\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   rc = IV_writeToBinaryFile(&inpmtx->sizesIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from writing sizes\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
   rc = IV_writeToBinaryFile(&inpmtx->offsetsIV, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in InpMtx_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from writing offsets\n", 
              inpmtx, fp, rc) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to write a InpMtx object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   ----------------------------------------------------
*/
int
InpMtx_writeForHumanEye ( 
   InpMtx    *inpmtx, 
   FILE   *fp 
) {
double   *entries ;
int      ierr, ii, iv, rc, size, vecid ;
int      *indices ;

if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_writeForHumanEye(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   exit(-1) ;
}
/*
   ------------------------
   write out the statistics
   ------------------------
*/
if ( (rc = InpMtx_writeStats(inpmtx, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in InpMtx_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from InpMtx_writeStats(%p,%p)\n",
           inpmtx, fp, rc, inpmtx, fp) ;
   return(0) ;
}
if (  inpmtx->nent > 0 ) {
   if ( INPMTX_IS_RAW_DATA(inpmtx) || INPMTX_IS_SORTED(inpmtx) ) {
      int   *ivec1 = InpMtx_ivec1(inpmtx) ;
      int   *ivec2 = InpMtx_ivec2(inpmtx) ;

      fprintf(fp, "\n data via triples") ;
      if ( INPMTX_IS_INDICES_ONLY(inpmtx) ) {
         for ( ii = 0 ; ii < inpmtx->nent ; ii++ ) {
            if ( ii % 4 == 0 ) {   
               fprintf(fp, "\n") ; 
            }
            fprintf(fp, " <%6d,%6d>", ivec1[ii], ivec2[ii]) ;
         }
      } else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
         double   *dvec = InpMtx_dvec(inpmtx) ;
         for ( ii = 0 ; ii < inpmtx->nent ; ii++ ) {
            if ( ii % 2 == 0 ) { fprintf(fp, "\n") ; }
            fprintf(fp, " <%6d,%6d,%20.12e>", 
                    ivec1[ii], ivec2[ii], dvec[ii]) ;
         }
      } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
         double   *dvec = InpMtx_dvec(inpmtx) ;
         for ( ii = 0 ; ii < inpmtx->nent ; ii++ ) {
            fprintf(fp, "\n <%6d,%6d,%20.12e,%20.12e>", 
                    ivec1[ii], ivec2[ii], dvec[2*ii], dvec[2*ii+1]) ;
         }
      }
   } else if ( INPMTX_IS_BY_VECTORS(inpmtx) && inpmtx->nvector > 0 ) {
      int   *vecids  = InpMtx_vecids(inpmtx) ;

      fprintf(fp, "\n data via vectors") ;
      if ( INPMTX_IS_INDICES_ONLY(inpmtx) ) {
         for ( iv = 0 ; iv < inpmtx->nvector ; iv++ ) {
            vecid = vecids[iv] ;
            InpMtx_vector(inpmtx, vecid, &size, &indices) ;
            if ( size > 0 ) {
               fprintf(fp, "\n %6d : ", vecids[iv]) ;
               IVfp80(fp, size, indices, 10, &ierr) ;
            }
         }
      } else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
         for ( iv = 0 ; iv < inpmtx->nvector ; iv++ ) {
            vecid = vecids[iv] ;
            InpMtx_realVector(inpmtx, vecid, &size, &indices, &entries);
            fprintf(fp, "\n %6d : ", vecids[iv]) ;
            IVfp80(fp, size, indices, 10, &ierr) ;
            DVfprintf(fp, size, entries) ;
         }
      } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
         for ( iv = 0 ; iv < inpmtx->nvector ; iv++ ) {
            vecid = vecids[iv] ;
            InpMtx_complexVector(inpmtx, vecid, 
                                 &size, &indices, &entries);
            fprintf(fp, "\n %6d : ", vecids[iv]) ;
            IVfp80(fp, size, indices, 10, &ierr) ;
            ZVfprintf(fp, size, entries) ;
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the InpMtx object

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   -----------------------------------------------------------
*/
int
InpMtx_writeStats ( 
   InpMtx    *inpmtx, 
   FILE   *fp 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in InpMtx_writeStats(%p,%p)"
           "\n bad input\n", inpmtx, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n InpMtx : double precision input matrix object :") ;
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   fprintf(fp, "\n coordType --> row triples") ;
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   fprintf(fp, "\n coordType --> column triples") ;
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   fprintf(fp, "\n coordType --> chevron triples") ;
} else if ( INPMTX_IS_CUSTOM(inpmtx) ) {
   fprintf(fp, "\n coordType --> custom triples") ;
} else {
   fprintf(stderr, "\n fatal error in InpMtx_writeStats(%p,%p)"
           "\n invalid inpmtx->coordType = %d\n", 
           inpmtx, fp, inpmtx->coordType) ;
   return(0) ;
}
if ( INPMTX_IS_RAW_DATA(inpmtx) ) {
   fprintf(fp, "\n storageMode --> raw triples") ;
} else if ( INPMTX_IS_SORTED(inpmtx) ) {
   fprintf(fp, "\n storageMode --> sorted and distinct triples") ;
} else if ( INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(fp, "\n storageMode --> vectors by first coordinate") ;
} else {
   fprintf(stderr, "\n fatal error in InpMtx_writeStats(%p,%p)"
           "\n invalid inpmtx->storageMode = %d\n", 
           inpmtx, fp, inpmtx->storageMode) ;
   return(0) ;
}
if ( INPMTX_IS_INDICES_ONLY(inpmtx) ) {
   fprintf(fp, "\n inputMode --> indices only") ;
} else if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
   fprintf(fp, "\n inputMode --> real entries") ;
} else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
   fprintf(fp, "\n inputMode --> complex entries") ;
} else {
   fprintf(stderr, "\n fatal error in InpMtx_writeStats(%p,%p)"
           "\n invalid inpmtx->inputMode = %d\n", 
           inpmtx, fp, inpmtx->inputMode) ;
   return(0) ;
}
fprintf(fp, "\n maxnent = %d --> maximum number of entries",
             inpmtx->maxnent) ;
fprintf(fp, "\n nent = %d --> present number of entries",
             inpmtx->nent) ;
fprintf(fp, "\n resizeMultiple = %.4g --> resize multiple",
             inpmtx->resizeMultiple) ;
fprintf(fp, "\n maxnvector = %d --> maximum number of vectors",
             inpmtx->maxnvector) ;
fprintf(fp, "\n nvector = %d --> present number of vectors",
             inpmtx->nvector) ;
fflush(fp) ;
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to write a InpMtx object to a matlab file

   return value -- 1 if success, 0 otherwise

   created -- 98jan28, cca
   ----------------------------------------------------
*/
int
InpMtx_writeForMatlab ( 
   InpMtx   *inpmtx, 
   char     *mtxname,
   FILE     *fp 
) {
int      ii, oldCoordType, oldStorageMode ;

if ( inpmtx == NULL || mtxname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in InpMtx_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", inpmtx, mtxname, fp) ;
   exit(-1) ;
}
/*
fprintf(fp, "\n before changing") ;
InpMtx_writeForHumanEye(inpmtx, fp) ;
*/
oldCoordType = inpmtx->coordType ;
oldStorageMode = inpmtx->storageMode ;
if ( oldCoordType != INPMTX_BY_ROWS ) {
   InpMtx_changeCoordType(inpmtx, INPMTX_BY_ROWS) ;
}
/*
fprintf(fp, "\n after changing") ;
InpMtx_writeForHumanEye(inpmtx, fp) ;
*/
if (  inpmtx->nent > 0 ) {
   int   *ivec1 = InpMtx_ivec1(inpmtx) ;
   int   *ivec2 = InpMtx_ivec2(inpmtx) ;

   if ( INPMTX_IS_REAL_ENTRIES(inpmtx) ) {
      double   *dvec = InpMtx_dvec(inpmtx) ;
      for ( ii = 0 ; ii < inpmtx->nent ; ii++ ) {
         fprintf(fp, "\n %s(%d,%d) = %24.16e ;", 
                 mtxname, 1+ivec1[ii], 1+ivec2[ii], dvec[ii]) ;
      }
   } else if ( INPMTX_IS_COMPLEX_ENTRIES(inpmtx) ) {
      double   *dvec = InpMtx_dvec(inpmtx) ;
      for ( ii = 0 ; ii < inpmtx->nent ; ii++ ) {
         fprintf(fp, "\n %s(%d,%d) = %24.16e + %24.16e*i;", 
                 mtxname, 1+ivec1[ii], 1+ivec2[ii], 
                 dvec[2*ii], dvec[2*ii+1]) ;
      }
   }
}
if ( oldCoordType != INPMTX_BY_ROWS ) {
   InpMtx_changeCoordType(inpmtx, oldCoordType) ;
}
InpMtx_changeStorageMode(inpmtx, oldStorageMode) ;
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- to read in a InpMtx object from a Harwell-Boeing file

   input --

      fn -- filename

   return value -- 1 if success, 0 if failure

   created -- 98sep11, cca
   ---------------------------------------------------------------
*/
int
InpMtx_readFromHBfile ( 
   InpMtx   *inpmtx, 
   char     *fn 
) {
char     *type ;
double   *values ;
int      ii, iiend, iistart, inputMode, jcol, ncol, nnonzeros,
         nrhs, nrow, rc ;
int      *colind, *colptr, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( inpmtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in InpMtx_readFromFile(%p,%s)"
           "\n bad input\n", inpmtx, fn) ;
   return(0) ;
}
/*
   ---------------------------------------------
   read in the Harwell-Boeing matrix information
   ---------------------------------------------
*/
if ( strcmp(fn, "none") == 0 ) {
   fprintf(stderr, "\n no file to read from") ;
   exit(0) ;
}
rc = readHB_info(fn, &nrow, &ncol, &nnonzeros, &type, &nrhs) ;
if ( rc != 1 ) {
   return(rc) ;
}
switch ( type[0] ) {
case 'P' :
   inputMode = INPMTX_INDICES_ONLY ;
   break ;
case 'R' :
   inputMode = SPOOLES_REAL ;
   break ;
case 'C' :
   inputMode = SPOOLES_COMPLEX ;
   break ;
default :
   fprintf(stderr, "\n fatal error in InpMtx_readFromHBfile"
           "\n type = %s, first character must be 'P', 'R' or 'C'",
           type) ;
   exit(-1) ;
   break ;
}
/*
   -----------------------------
   initialize the InpMtx object
   -----------------------------
*/
InpMtx_init(inpmtx, INPMTX_BY_COLUMNS, inputMode, nnonzeros, 0) ;
colptr = IVinit(ncol+1, -1) ;
colind = InpMtx_ivec1(inpmtx)   ;
rowind = InpMtx_ivec2(inpmtx)   ;
values = InpMtx_dvec(inpmtx)    ;
InpMtx_setNent(inpmtx, nnonzeros) ;
/*
   -------------------------------
   read in the indices and entries
   -------------------------------
*/
rc = readHB_mat_double(fn, colptr, rowind, values) ;
if ( rc != 1 ) {
   return(rc) ;
}
/*
   --------------------------------------------
   decrement the column offsets and row indices
   --------------------------------------------
*/
for ( jcol = 0 ; jcol <= ncol ; jcol++ ) {
   colptr[jcol]-- ;
}
for ( ii = 0 ; ii < nnonzeros ; ii++ ) {
   rowind[ii]-- ;
}
/*
   -------------------------------------------
   fill the ivec1[] vector with column indices
   -------------------------------------------
*/
for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
   iistart = colptr[jcol] ;
   iiend   = colptr[jcol+1] - 1 ;
   for ( ii = iistart ; ii <= iiend ; ii++ ) {
      colind[ii] = jcol ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(colptr) ;
CVfree(type) ;

return(1) ; }

/*--------------------------------------------------------------------*/
