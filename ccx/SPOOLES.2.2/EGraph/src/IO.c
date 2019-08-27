/*  IO.c  */

#include "../EGraph.h"

static const char *suffixb = ".egraphb" ;
static const char *suffixf = ".egraphf" ;

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to read in a EGraph object from a file

   input --

      fn -- filename, must be *.egraphb or *.egraphf

   return value -- 1 if success, 0 if failure

   created -- 95nov03, cca
   -------------------------------------------------
*/
int
EGraph_readFromFile ( 
   EGraph   *egraph, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in EGraph_readFromFile(%p,%s)"
           "\n bad input\n", egraph, fn) ;
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
         fprintf(stderr, "\n error in EGraph_readFromFile(%p,%s)"
                 "\n unable to open file %s", egraph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = EGraph_readFromBinaryFile(egraph, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in EGraph_readFromFile(%p,%s)"
                 "\n unable to open file %s", egraph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = EGraph_readFromFormattedFile(egraph, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in EGraph_readFromFile(%p,%s)"
              "\n bad EGraph file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              egraph, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in EGraph_readFromFile(%p,%s)"
       "\n bad EGraph file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       egraph, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read a EGraph object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95nov03, cca
   --------------------------------------------------------
*/
int
EGraph_readFromFormattedFile ( 
   EGraph   *egraph, 
   FILE     *fp 
) {
int   nelem, nvtx, rc, type ;
int   itemp[6] ;
int   *vwghts ;
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in EGraph_readFromFormattedFile(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
EGraph_clearData(egraph) ;
/*
   -----------------------------------
   read in the three scalar parameters
   type, nelem, nvtx
   -----------------------------------
*/
if ( (rc = IVfscanf(fp, 3, itemp)) != 3 ) {
   fprintf(stderr, "\n error in EGraph_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", egraph, fp, rc, 3) ;
   return(0) ;
}
type  = itemp[0] ;
nelem = itemp[1] ;
nvtx  = itemp[2] ;
/*
   ----------------------------
   initialize the EGraph object
   ----------------------------
*/
EGraph_init(egraph, type, nelem, nvtx, IVL_CHUNKED) ;
/*
   --------------------------
   read in the adj IVL object
   --------------------------
*/
adjIVL = egraph->adjIVL ;
rc = IVL_readFromFormattedFile(adjIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in EGraph_readFromFormattedFile(%p,%p)"
           "\n trying to read in adjIVL"
           "\n return code %d from IVL_readFormattedFile(%p,%p)",
           egraph, fp, rc, adjIVL, fp) ;
   return(0) ;
}
if ( type % 2 == 1 ) {
/*
   --------------------------
   vertex weights are present
   --------------------------
*/
   vwghts = egraph->vwghts ;
   if ( (rc = IVfscanf(fp, nvtx, vwghts)) != nvtx ) {
      fprintf(stderr, "\n error in EGraph_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", egraph, fp, rc, nvtx) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to read a EGraph object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95nov03, cca
   ----------------------------------------------------
*/
int
EGraph_readFromBinaryFile ( 
   EGraph   *egraph, 
   FILE    *fp 
) {
int   nelem, nvtx, rc, type ;
int   itemp[6] ;
int   *vwghts ;
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in EGraph_readFromBinaryFile(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
EGraph_clearData(egraph) ;
/*
   -----------------------------------
   read in the three scalar parameters
   type, nelem, nvtx
   -----------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 3, fp)) != 3 ) {
   fprintf(stderr, "\n error in EGraph_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", egraph, fp, rc, 3) ;
   return(0) ;
}
type  = itemp[0] ;
nelem = itemp[1] ;
nvtx  = itemp[2] ;
/*
   ----------------------------
   initialize the EGraph object
   ----------------------------
*/
EGraph_init(egraph, type, nelem, nvtx, IVL_CHUNKED) ;
/*
   --------------------------
   read in the adj IVL object
   --------------------------
*/
adjIVL = egraph->adjIVL ;
rc = IVL_readFromBinaryFile(adjIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in EGraph_readFromBinaryFile(%p,%p)"
           "\n trying to read in adjIVL"
           "\n return code %d from IVL_readBinaryFile(%p,%p)",
           egraph, fp, rc, adjIVL, fp) ;
   return(0) ;
}
if ( type % 2 == 1 ) {
/*
   --------------------------
   vertex weights are present
   --------------------------
*/
   vwghts = egraph->vwghts ;
   if ( (rc = fread((void *) vwghts, sizeof(int), nvtx, fp)) != nvtx){
      fprintf(stderr, "\n error in EGraph_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", egraph, fp, rc, nvtx) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to write a EGraph object to a file

   input --

      fn -- filename
        *.egraphb -- binary
        *.egraphf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov03, cca
   --------------------------------------------
*/
int
EGraph_writeToFile ( 
   EGraph   *egraph, 
   char     *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in EGraph_writeToFile(%p,%s)"
    "\n bad input\n", egraph, fn) ; 
   return(0) ;
}
if ( egraph->type < 0 || 3 < egraph->type ) {
   fprintf(stderr, "\n fatal error in EGraph_writeToFile(%p,%s)"
           "\n bad type = %d", egraph, fn, egraph->type) ;
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
         fprintf(stderr, "\n error in EGraph_writeToFile(%p,%s)"
                 "\n unable to open file %s", egraph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = EGraph_writeToBinaryFile(egraph, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in EGraph_writeToFile(%p,%s)"
                 "\n unable to open file %s", egraph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = EGraph_writeToFormattedFile(egraph, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in EGraph_writeToFile(%p,%s)"
                 "\n unable to open file %s", egraph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = EGraph_writeForHumanEye(egraph, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in EGraph_writeToFile(%p,%s)"
              "\n unable to open file %s", egraph, fn, fn) ;
      rc = 0 ;
   } else {
      rc = EGraph_writeForHumanEye(egraph, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write a EGraph object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95nov03, cca
   ------------------------------------------------------
*/
int
EGraph_writeToFormattedFile ( 
   EGraph   *egraph, 
   FILE     *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   return(0) ;
}
if ( egraph->type < 0 || 1 < egraph->type ) {
   fprintf(stderr, 
           "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
           "\n bad type = %d", egraph, fp, egraph->type) ;
   return(0) ;
}
/*
   -------------------------------------
   write out the three scalar parameters
   -------------------------------------
*/
rc = fprintf(fp, "\n %d %d %d", 
             egraph->type, egraph->nelem, egraph->nvtx) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", egraph, fp, rc) ;
   return(0) ;
}
/*
   ---------------------------------
   write out the adjacency structure
   ---------------------------------
*/
rc = IVL_writeToFormattedFile(egraph->adjIVL, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from IVL_writeToFormattedFile(%p,%p)" 
           "\n while attempting to write out adjIVL\n",
           egraph, fp, rc, egraph->adjIVL, fp) ;
   return(0) ;
}
/*
   -----------------------------------------
   write out the vwghts[] vector, if present
   -----------------------------------------
*/
if ( egraph->type % 2 == 1 ) {
   if ( egraph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
              "\n egraph->type = %d, egraph->vwghts == NULL\n",
              egraph, fp, egraph->type) ;
      return(0) ;
   }
   IVfp80(fp, egraph->nvtx, egraph->vwghts, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from vwghts[] IVfp80\n", 
              egraph, fp, ierr) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write a EGraph object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95nov03, cca
   ---------------------------------------------------
*/
int
EGraph_writeToBinaryFile ( 
   EGraph   *egraph, 
   FILE     *fp 
) {
int   rc ;
int   itemp[6] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in EGraph_writeToBinaryFile(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   return(0) ;
}
if ( egraph->type < 0 || 3 < egraph->type ) {
   fprintf(stderr, "\n fatal error in EGraph_writeToBinaryFile(%p,%p)"
           "\n bad type = %d", egraph, fp, egraph->type) ;
   return(0) ;
}
/*
   -------------------------------------
   write out the three scalar parameters
   -------------------------------------
*/
itemp[0] = egraph->type  ;
itemp[1] = egraph->nelem ;
itemp[2] = egraph->nvtx  ;
rc = fwrite((void *) itemp, sizeof(int), 6, fp) ;
if ( rc != 6 ) {
   fprintf(stderr, "\n error in EGraph_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", egraph, fp, rc, 6) ;
   return(0) ;
}
/*
   ---------------------------------
   write out the adjacency structure
   ---------------------------------
*/
rc = IVL_writeToBinaryFile(egraph->adjIVL, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in EGraph_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from IVL_writeToBinaryFile(%p,%p)" 
           "\n while attempting to write out adj\n",
           egraph, fp, rc, egraph->adjIVL, fp) ;
   return(0) ;
}
/*
   -----------------------------------------
   write out the vwghts[] vector, if present
   -----------------------------------------
*/
if ( egraph->type % 2 == 1 ) {
   if ( egraph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeToBinaryFile(%p,%p)"
              "\n egraph->type = %d, egraph->vwghts == NULL\n",
              egraph, fp, egraph->type) ;
      return(0) ;
   }
   rc = fwrite((void *) egraph->vwghts, sizeof(int), egraph->nvtx, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from vwghts[] fwrite\n", 
              egraph, fp, rc) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write a EGraph object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov03, cca
   -------------------------------------------------
*/
int
EGraph_writeForHumanEye ( 
   EGraph    *egraph, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in EGraph_writeForHumanEye(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   exit(-1) ;
}
/*
   ------------------------
   write out the statistics
   ------------------------
*/
if ( (rc = EGraph_writeStats(egraph, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in EGraph_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from EGraph_writeStats(%p,%p)\n",
           egraph, fp, rc, egraph, fp) ;
   return(0) ;
}
if ( egraph->adjIVL != NULL ) {
/*
   ----------------------------------
   write out the adjacency IVL object
   ----------------------------------
*/
   fprintf(fp, "\n\n adjacency IVL object") ;
   rc = IVL_writeForHumanEye(egraph->adjIVL, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, "\n fatal error in EGraph_writeForHumanEye(%p,%p)"
              "\n rc = %d, return from IVL_writeForHumanEye(%p,%p)" 
              "\n while attempting to write out adjIVL\n",
              egraph, fp, rc, egraph->adjIVL, fp) ;
      return(0) ;
   }
}
/*
   ----------------------------------------------
   write out the vertex weights vector if present
   ----------------------------------------------
*/
if ( egraph->type % 2 == 1 ) {
   if ( egraph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeForHumanEye(%p,%p)"
              "\n egraph->type = %d, egraph->vwghts == NULL\n",
              egraph, fp, egraph->type) ;
      return(0) ;
   }
   fprintf(fp, "\n\n vertex weights ") ;
   IVfp80(fp, egraph->nvtx, egraph->vwghts, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in EGraph_writeForHumanEye(%p,%p)"
              "\n ierr = %d, return from vwghts[] IVfp80\n", 
              egraph, fp, ierr) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the EGraph object

   return value -- 1 if success, 0 otherwise

   created -- 95nov03, cca
   -----------------------------------------------------------
*/
int
EGraph_writeStats ( 
   EGraph    *egraph, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( egraph == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in EGraph_writeStats(%p,%p)"
           "\n bad input\n", egraph, fp) ;
   exit(-1) ;
}
switch ( egraph->type ) {
case 0 : 
   rc = fprintf(fp, "\n EGraph : unweighted egraph object :") ;
   break ;
case 1 : 
   rc = fprintf(fp, "\n EGraph : vertices weighted egraph object :") ;
   break ;
default :
   fprintf(stderr, "\n fatal error in EGraph_writeStats(%p,%p)"
           "\n invalid egraph->type = %d\n", egraph, fp, egraph->type) ;
   return(0) ;
}
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;
rc = fprintf(fp, " %d elements, %d vertices",
             egraph->nelem, egraph->nvtx) ;
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;

return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in EGraph_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", egraph, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
