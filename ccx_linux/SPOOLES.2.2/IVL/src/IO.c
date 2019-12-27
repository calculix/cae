/*  IO.c  */

#include "../IVL.h"

static const char *suffixb = ".ivlb" ;
static const char *suffixf = ".ivlf" ;

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to read in an IVL object from a file

   input --

      fn -- filename, must be *.ivlb or *.ivlf

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   -----------------------------------------------
*/
int
IVL_readFromFile ( 
   IVL    *ivl, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in IVL_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", ivl, fn, __FILE__, __LINE__) ;
   return(0) ;
}
switch ( ivl->type ) {
case IVL_CHUNKED :
case IVL_SOLO    :
case IVL_UNKNOWN :
   break ;
default :
   fprintf(stderr, "\n error in IVL_readFromFile(%p,%s)"
    "\n bad type = %d", ivl, fn, ivl->type) ;
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
         fprintf(stderr, "\n error in IVL_readFromFile(%p,%s)"
                 "\n unable to open file %s", ivl, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IVL_readFromBinaryFile(ivl, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in IVL_readFromFile(%p,%s)"
                 "\n unable to open file %s", ivl, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IVL_readFromFormattedFile(ivl, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in IVL_readFromFile(%p,%s)"
              "\n bad IVL file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              ivl, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in IVL_readFromFile(%p,%s)"
       "\n bad IVL file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       ivl, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read an IVL object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   ------------------------------------------------------
*/
int
IVL_readFromFormattedFile ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   nlist, rc, type ;
int   itemp[3] ;
int   *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in IVL_readFromFormattedFile(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   return(0) ;
}
switch ( ivl->type ) {
case IVL_CHUNKED :
case IVL_SOLO    :
   break ;
default :
   fprintf(stderr, "\n error in IVL_readFormattedFile(%p,%p)"
    "\n bad type = %d", ivl, fp, ivl->type) ;
   return(0) ;
}
/*
   -------------------------------------------
   save the ivl type and clear the data fields
   -------------------------------------------
*/
type = ivl->type ;
IVL_clearData(ivl) ;
/*
   -----------------------------------
   read in the three scalar parameters
   type, # of lists, # of indices
   -----------------------------------
*/
if ( (rc = IVfscanf(fp, 3, itemp)) != 3 ) {
   fprintf(stderr, "\n error in IVL_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", ivl, fp, rc, 3) ;
   return(0) ;
}
nlist = itemp[1] ;
/*
fprintf(stdout, "\n itemp = { %d %d %d } ",
        itemp[0], itemp[1], itemp[2]) ;
*/
/*
   --------------------------
   read in the sizes[] vector
   --------------------------
*/
sizes = IVinit(nlist, 0) ;
if ( (rc = IVfscanf(fp, nlist, sizes)) != nlist ) {
   fprintf(stderr, "\n error in IVL_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", ivl, fp, rc, nlist) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
IVL_init3(ivl, type, nlist, sizes) ;
IVfree(sizes) ;
/*
   -----------------------
   now read in the indices
   -----------------------
*/
switch ( type ) {
case IVL_SOLO : {
   int   ilist, size ;
   int   *ind ;

   for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
      IVL_listAndSize(ivl, ilist, &size, &ind) ;
      if ( size > 0 ) {
         if ( (rc = IVfscanf(fp, size, ind)) != size ) {
            fprintf(stderr, 
                    "\n error in IVL_readFromFormattedFile(%p,%p)"
                    "\n list %d, %d items of %d read\n", 
                    ivl, fp, ilist, rc, size) ;
            return(0) ;
         }
      }
   }
   } break ;
case IVL_CHUNKED : {
/*
   --------------------------------------------------------
   read in the indices into the contiguous block of storage
   --------------------------------------------------------
*/
   if ( (rc = IVfscanf(fp, ivl->tsize, ivl->chunk->base)) 
        != ivl->tsize ) {
      fprintf(stderr, "\n error in IVL_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", ivl, fp, rc, ivl->tsize) ;
      return(0) ;
   }
   } break ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to read an IVL object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95sep29, cca
   ---------------------------------------------------
*/
int
IVL_readFromBinaryFile ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   nlist, rc, type ;
int   itemp[3] ;
int   *sizes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_readFromBinaryFile(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   return(0) ;
}
switch ( ivl->type ) {
case IVL_CHUNKED :
case IVL_SOLO    :
   break ;
default :
   fprintf(stderr, "\n error in IVL_readBinaryFile(%p,%p)"
    "\n bad type = %d", ivl, fp, ivl->type) ;
   return(0) ;
}
/*
   -------------------------------------------
   save the ivl type and clear the data fields
   -------------------------------------------
*/
type = ivl->type ;
IVL_clearData(ivl) ;
/*
   -----------------------------------
   read in the three scalar parameters
   type, # of lists, # of indices
   -----------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 3, fp)) != 3 ) {
   fprintf(stderr, "\n error in IVL_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", ivl, fp, rc, 3) ;
   return(0) ;
}
nlist = itemp[1] ;
/*
   --------------------------
   read in the sizes[] vector
   --------------------------
*/
sizes = IVinit(nlist, 0) ;
if ( (rc = fread((void *) sizes, sizeof(int), nlist, fp)) != nlist ) {
   fprintf(stderr, "\n error in IVL_readFromBinaryFile(%p,%p)"
           "\n sizes(%d) : %d items of %d read\n", 
           ivl, fp, nlist, rc, nlist) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
IVL_init3(ivl, type, nlist, sizes) ;
IVfree(sizes) ;
/*
   -------------------
   read in the indices
   -------------------
*/
switch ( type ) {
case IVL_SOLO : {
   int   ilist, size ;
   int   *ind ;

   for ( ilist = 0 ; ilist < nlist ; ilist++ ) {
      IVL_listAndSize(ivl, ilist, &size, &ind) ;
      if ( (rc = fread((void *) ind, sizeof(int), size, fp)) != size ) {
         fprintf(stderr, "\n error in IVL_readFromBinaryFile(%p,%p)"
                 "\n list %d, %d items of %d read\n", 
                 ivl, fp, ilist, rc, size) ;
         return(0) ;
      }
   }
   } break ;
case IVL_CHUNKED : {
/*
   --------------------------------------------------------
   read in the indices into the contiguous block of storage
   --------------------------------------------------------
*/
   if ( (rc = fread((void *) ivl->chunk->base, 
                    sizeof(int), ivl->tsize, fp)) != ivl->tsize ) {
      fprintf(stderr, "\n error in IVL_readFromBinaryFile(%p,%p)"
              "\n indices(%d) : %d items of %d read\n", 
              ivl, fp, ivl->tsize, rc, ivl->tsize) ;
      return(0) ;
   }
   } break ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   purpose -- to write an IVL object to a file

   input --

      fn -- filename
        *.ivlb -- binary
        *.ivlf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -------------------------------------------
*/
int
IVL_writeToFile ( 
   IVL    *ivl, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_writeToFile(%p,%s)"
    "\n bad input\n", ivl, fn) ; 
}
switch ( ivl->type ) {
case IVL_CHUNKED :
case IVL_SOLO    :
case IVL_UNKNOWN :
   break ;
default :
   fprintf(stderr, "\n fatal error in IVL_writeToFile(%p,%s)"
           "\n bad type = %d", ivl, fn, ivl->type) ;
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
         fprintf(stderr, "\n error in IVL_writeToFile(%p,%s)"
                 "\n unable to open file %s", ivl, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IVL_writeToBinaryFile(ivl, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in IVL_writeToFile(%p,%s)"
                 "\n unable to open file %s", ivl, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IVL_writeToFormattedFile(ivl, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in IVL_writeToFile(%p,%s)"
                 "\n unable to open file %s", ivl, fn, fn) ;
         rc = 0 ;
      } else {
         rc = IVL_writeForHumanEye(ivl, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in IVL_writeToFile(%p,%s)"
              "\n unable to open file %s", ivl, fn, fn) ;
      rc = 0 ;
   } else {
      rc = IVL_writeForHumanEye(ivl, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an IVL object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -----------------------------------------------------
*/
int
IVL_writeToFormattedFile ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   count, ierr, j, jsize, nlist, rc ;
int   *jind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fp == NULL || (nlist = ivl->nlist) <= 0 ) {
   fprintf(stderr, "\n fatal error in IVL_writeToFormattedFile(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the three scalar parameters
   -------------------------------------
*/
rc = fprintf(fp, "\n %d %d %d", ivl->type, ivl->nlist, ivl->tsize) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in IVL_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", ivl, fp, rc) ;
   return(0) ;
}
if ( ivl->nlist > 0 ) {
   IVfp80(fp, ivl->nlist, ivl->sizes, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in IVL_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from sizes[] IVfp80\n", 
              ivl, fp, ierr) ;
      return(0) ;
   }
}
switch ( ivl->type ) {
case IVL_NOTYPE :
   break ;
case IVL_UNKNOWN :
case IVL_CHUNKED :
case IVL_SOLO    :
   for ( j = 0, count = 80 ; j < ivl->nlist ; j++ ) {
      IVL_listAndSize(ivl, j, &jsize, &jind) ;
      if ( jsize > 0 ) {
         count = IVfp80(fp, jsize, jind, count, &ierr) ;
         if ( ierr < 0 ) {
            fprintf(stderr, 
                    "\n fatal error in IVL_writeToFormattedFile(%p,%p)"
                    "\n ierr = %d, return from IVfp80, list %d\n", 
                    ivl, fp, ierr, j) ;
            return(0) ;
         }
      }
   }
   break ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an IVL object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   --------------------------------------------------
*/
int
IVL_writeToBinaryFile ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   j, jsize, nlist, rc ;
int   *jind ;
int   itemp[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fp == NULL || (nlist = ivl->nlist) <= 0 ) {
   fprintf(stderr, "\n fatal error in IVL_writeToBinaryFile(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   exit(-1) ;
}
itemp[0] = ivl->type ;
itemp[1] = ivl->nlist ;
itemp[2] = ivl->tsize ;
rc = fwrite((void *) itemp, sizeof(int), 3, fp) ;
if ( rc != 3 ) {
   fprintf(stderr, "\n error in IVL_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", ivl, fp, rc, 3) ;
   return(0) ;
}
rc = fwrite((void *) ivl->sizes, sizeof(int), ivl->nlist, fp) ;
if ( rc != ivl->nlist ) {
   fprintf(stderr, "\n error in IVL_writeToBinaryFile(%p,%p)"
           "\n ivl->sizes, %d of %d items written\n",
           ivl, fp, rc, ivl->nlist) ;
   return(0) ;
}
switch ( ivl->type ) {
case IVL_NOTYPE :
   break ;
case IVL_CHUNKED :
case IVL_SOLO    :
case IVL_UNKNOWN :
   for ( j = 0 ; j < ivl->nlist ; j++ ) {
      IVL_listAndSize(ivl, j, &jsize, &jind) ;
      if ( jsize > 0 ) {
         rc = fwrite((void *) jind, sizeof(int), jsize, fp) ;
         if ( rc != jsize ) {
            fprintf(stderr, "\n error in IVL_writeToBinaryFile(%p,%p)"
                    "\n list %d, %d of %d items written\n",
                    ivl, fp, j, rc, jsize) ;
            return(0) ;
         }
      }
   }
   break ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an IVL object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -------------------------------------------------
*/
int
IVL_writeForHumanEye ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   ierr, j, size, rc ;
int   *ind ;

if ( ivl == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in IVL_writeForHumanEye(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   exit(-1) ;
}
if ( (rc = IVL_writeStats(ivl, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in IVL_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from IVL_writeStats(%p,%p)\n",
           ivl, fp, rc, ivl, fp) ;
   return(0) ;
}
for ( j = 0 ; j < ivl->nlist ; j++ ) {
   IVL_listAndSize(ivl, j, &size, &ind) ;
   if ( size > 0 ) {
      fprintf(fp, "\n %5d :", j) ;
      IVfp80(fp, size, ind, 8, &ierr) ;
      if ( ierr < 0 ) {
         fprintf(stderr, "\n fatal error in IVL_writeForHumanEye(%p,%p)"
                 "\n ierr = %d, return from IVfp80, list %d\n", 
                 ivl, fp, ierr, j) ;
         return(0) ;
      }
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to write out the statistics for the IVL object

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   ---------------------------------------------------------
*/
int
IVL_writeStats ( 
   IVL    *ivl, 
   FILE   *fp 
) {
int   nactive, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( ivl == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in IVL_writeStats(%p,%p)"
           "\n bad input\n", ivl, fp) ;
   exit(-1) ;
}
nactive = 0 ;
if ( ivl->nlist > 0 ) {
   nactive = IVsum(ivl->nlist, ivl->sizes) ;
}
rc = fprintf(fp, "\n IVL : integer vector list object :") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, "\n type %d", ivl->type) ;
if ( rc < 0 ) { goto IO_error ; }
switch ( ivl->type ) {
case IVL_CHUNKED : rc = fprintf(fp, ", chunked storage") ; break ;
case IVL_SOLO    : rc = fprintf(fp, ", solo storage")    ; break ;
case IVL_UNKNOWN : rc = fprintf(fp, ", unknown storage") ; break ;
}
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, 
             "\n %d lists, %d maximum lists, %d tsize, %d total bytes",
             ivl->nlist, ivl->maxnlist, ivl->tsize, IVL_sizeOf(ivl)) ;
if ( rc < 0 ) { goto IO_error ; }
switch ( ivl->type ) {
case IVL_CHUNKED : {
   Ichunk   *chunk ;
   int      nalloc, nchunk ;

   nalloc = nchunk = 0 ;
   for ( chunk = ivl->chunk ; chunk != NULL ; chunk = chunk->next){
      nchunk++ ;
      nalloc += chunk->size ;
   }
   rc = fprintf(fp, "\n %d chunks, %d active entries, %d allocated",
                nchunk, nactive, nalloc) ;
   if ( rc < 0 ) { goto IO_error ; }
   if ( nalloc > 0 ) {
      rc = fprintf(fp, ", %.2f %% used", (100.*nactive)/nalloc) ;
      if ( rc < 0 ) { goto IO_error ; }
   }
   } break ;
case IVL_SOLO :
   rc = fprintf(fp, 
                "\n %d lists separately allocated, %d active entries",
                ivl->nlist, nactive) ;
  if ( rc < 0 ) { goto IO_error ; }
   break ;
}
return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in IVL_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", ivl, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
