/*  IO.c  */

#include "../Coords.h"

static const char *suffixb = ".coordsb" ;
static const char *suffixf = ".coordsf" ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to read in an Coords object from a file

   input --

      fn -- filename, must be *.coordsb or *.coordsf

   return value -- 1 if success, 0 if failure

   created -- 95dec19, cca
   ---------------------------------------------------
*/
int
Coords_readFromFile ( 
   Coords   *coords, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || fn == NULL ) {
   fprintf(stderr, 
    "\n error in Coords_readFromFile(%p,%s), file %s, line %d"
    "\n bad input\n", coords, fn, __FILE__, __LINE__) ;
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
         fprintf(stderr, "\n error in Coords_readFromFile(%p,%s)"
                 "\n unable to open file %s", coords, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Coords_readFromBinaryFile(coords, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in Coords_readFromFile(%p,%s)"
                 "\n unable to open file %s", coords, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Coords_readFromFormattedFile(coords, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in Coords_readFromFile(%p,%s)"
              "\n bad Coords file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              coords, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in Coords_readFromFile(%p,%s)"
       "\n bad Coords file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       coords, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to read an Coords object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95dec19, cca
   ---------------------------------------------------------
*/
int
Coords_readFromFormattedFile ( 
   Coords   *coords, 
   FILE     *fp 
) {
int   rc, size ;
int   itemp[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Coords_readFromFormattedFile(%p,%p)"
           "\n bad input\n", coords, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Coords_clearData(coords) ;
/*
   ---------------------------------------
   read in the three scalar parameters
   type, # of dimensions, # of coordinates
   ---------------------------------------
*/
if ( (rc = IVfscanf(fp, 3, itemp)) != 3 ) {
   fprintf(stderr, "\n error in Coords_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", coords, fp, rc, 3) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Coords_init(coords, itemp[0], itemp[1], itemp[2]) ;
/*
   --------------------------
   read in the coors[] vector
   --------------------------
*/
size = itemp[1]*itemp[2] ;
if ( (rc = FVfscanf(fp, size, coords->coors)) != size ) {
   fprintf(stderr, "\n error in Coords_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", coords, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read an Coords object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95dec19, cca
   ------------------------------------------------------
*/
int
Coords_readFromBinaryFile ( 
   Coords   *coords, 
   FILE     *fp 
) {
int   rc, size ;
int   itemp[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_readFromBinaryFile(%p,%p)"
           "\n bad input\n", coords, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Coords_clearData(coords) ;
/*
   ---------------------------------------
   read in the three scalar parameters
   type, # of dimensions, # of coordinates
   ---------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 3, fp)) != 3 ) {
   fprintf(stderr, "\n error in Coords_readFromBinaryFile(%p,%p)"
           "\n itemp(3) : %d items of %d read\n", coords, fp, rc, 3) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Coords_init(coords, itemp[0], itemp[1], itemp[2]) ;
/*
   --------------------------
   read in the coors[] vector
   --------------------------
*/
size = itemp[1]*itemp[2] ;
if ( (rc = fread((void *) coords->coors, sizeof(float), size, fp)) 
      != size ) {
   fprintf(stderr, "\n error in Coords_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", coords, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to write an Coords object to a file

   input --

      fn -- filename
        *.coordsb -- binary
        *.coordsf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ----------------------------------------------
*/
int
Coords_writeToFile ( 
   Coords   *coords, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_writeToFile(%p,%s)"
    "\n bad input\n", coords, fn) ; 
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
         fprintf(stderr, "\n error in Coords_writeToFile(%p,%s)"
                 "\n unable to open file %s", coords, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Coords_writeToBinaryFile(coords, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in Coords_writeToFile(%p,%s)"
                 "\n unable to open file %s", coords, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Coords_writeToFormattedFile(coords, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in Coords_writeToFile(%p,%s)"
                 "\n unable to open file %s", coords, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Coords_writeForHumanEye(coords, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in Coords_writeToFile(%p,%s)"
              "\n unable to open file %s", coords, fn, fn) ;
      rc = 0 ;
   } else {
      rc = Coords_writeForHumanEye(coords, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to write an Coords object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   --------------------------------------------------------
*/
int
Coords_writeToFormattedFile ( 
   Coords   *coords, 
   FILE     *fp 
) {
int   rc, size ;
/*
   ---------------
   check the input
   ---------------
*/
if (   coords == NULL || fp == NULL 
   || (size = coords->ndim * coords->ncoor) <= 0 ) {
   fprintf(stderr, "\n fatal error in Coords_writeToFormattedFile(%p,%p)"
           "\n bad input\n", coords, fp) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the three scalar parameters
   -------------------------------------
*/
rc = fprintf(fp, "\n %d %d %d", 
             coords->type, coords->ndim, coords->ncoor) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in Coords_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", coords, fp, rc) ;
   return(0) ;
}
/*
   --------------------------------
   write out the coordinates vector
   --------------------------------
*/
FVfprintf(fp, size, coords->coors) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an Coords object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   -----------------------------------------------------
*/
int
Coords_writeToBinaryFile ( 
   Coords   *coords, 
   FILE     *fp 
) {
int   rc, size ;
int   itemp[3] ;
/*
   ---------------
   check the input
   ---------------
*/
if (   coords == NULL || fp == NULL 
   || (size = coords->ndim * coords->ncoor) <= 0 ) {
   fprintf(stderr, "\n fatal error in Coords_writeToBinaryFile(%p,%p)"
           "\n bad input\n", coords, fp) ;
   exit(-1) ;
}
itemp[0] = coords->type  ;
itemp[1] = coords->ndim  ;
itemp[2] = coords->ncoor ;
rc = fwrite((void *) itemp, sizeof(int), 3, fp) ;
if ( rc != 3 ) {
   fprintf(stderr, "\n error in Coords_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", coords, fp, rc, 3) ;
   return(0) ;
}
rc = fwrite((void *) coords->coors, sizeof(float), size, fp) ;
if ( rc != size ) {
   fprintf(stderr, "\n error in Coords_writeToBinaryFile(%p,%p)"
           "\n coords->coors, %d of %d items written\n",
           coords, fp, rc, size) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to write an Coords object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ----------------------------------------------------
*/
int
Coords_writeForHumanEye ( 
   Coords    *coords, 
   FILE   *fp 
) {
int   icoor, idim, rc ;

if ( coords == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Coords_writeForHumanEye(%p,%p)"
           "\n bad input\n", coords, fp) ;
   exit(-1) ;
}
if ( (rc = Coords_writeStats(coords, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in Coords_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from Coords_writeStats(%p,%p)\n",
           coords, fp, rc, coords, fp) ;
   return(0) ;
}
for ( icoor = 0 ; icoor < coords->ncoor ; icoor++ ) {
   fprintf(fp, "\n %6d", icoor) ;
   for ( idim = 1 ; idim <= coords->ndim ; idim++ ) {
      fprintf(fp, " %12.4g", Coords_value(coords, idim, icoor)) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to write out the statistics for the Coords object

   return value -- 1 if success, 0 otherwise

   created -- 95dec19, cca
   ------------------------------------------------------------
*/
int
Coords_writeStats ( 
   Coords    *coords, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( coords == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Coords_writeStats(%p,%p)"
           "\n bad input\n", coords, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n Coords : coordinates object :") ;
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, "\n          type %d", coords->type) ;
if ( rc < 0 ) { goto IO_error ; }
switch ( coords->type ) {
case 1 : rc = fprintf(fp, ", storage by tuples") ; break ;
case 2 : rc = fprintf(fp, ", storage by vectors") ; break ;
}
if ( rc < 0 ) { goto IO_error ; }
rc = fprintf(fp, 
      "\n          %d dimensions, %d coordinates, occupies %d bytes", 
      coords->ndim, coords->ncoor, Coords_sizeOf(coords)) ;

if ( rc < 0 ) { goto IO_error ; }

return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in Coords_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", coords, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
