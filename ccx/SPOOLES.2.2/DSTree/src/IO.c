/*  IO.c  */

#include "../DSTree.h"

static const char *suffixb = ".dstreeb" ;
static const char *suffixf = ".dstreef" ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to read in an DSTree object from a file

   input --

      fn -- filename, must be *.dstreeb or *.dstreef

   return value -- 1 if success, 0 if failure

   created -- 96mar10, cca
   --------------------------------------------------
*/
int
DSTree_readFromFile ( 
   DSTree   *dstree, 
   char    *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in DSTree_readFromFile(%p,%s)"
           "\n bad input\n", dstree, fn) ;
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
         fprintf(stderr, "\n error in DSTree_readFromFile(%p,%s)"
                 "\n unable to open file %s", dstree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DSTree_readFromBinaryFile(dstree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in DSTree_readFromFile(%p,%s)"
                 "\n unable to open file %s", dstree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DSTree_readFromFormattedFile(dstree, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in DSTree_readFromFile(%p,%s)"
              "\n bad DSTree file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              dstree, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in DSTree_readFromFile(%p,%s)"
       "\n bad DSTree file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       dstree, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- to read an DSTree object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 96mar10, cca
   ---------------------------------------------------------
*/
int
DSTree_readFromFormattedFile ( 
   DSTree   *dstree, 
   FILE    *fp 
) {
IV     *mapIV ;
Tree   *tree  ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DSTree_readFromFormattedFile(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
DSTree_clearData(dstree) ;
/*
   -----------------------
   read in the Tree object
   -----------------------
*/
tree = Tree_new() ;
Tree_readFromFormattedFile(tree, fp) ;
/*
   ---------------------
   read in the IV object
   ---------------------
*/
mapIV = IV_new() ;
IV_readFromFormattedFile(mapIV, fp) ;
/*
   ----------------------------
   initialize the DSTree object
   ----------------------------
*/
DSTree_init2(dstree, tree, mapIV) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to read an DSTree object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 96mar10, cca
   ------------------------------------------------------
*/
int
DSTree_readFromBinaryFile ( 
   DSTree    *dstree, 
   FILE   *fp 
) {
IV     *mapIV ;
Tree   *tree  ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_readFromBinaryFile(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
DSTree_clearData(dstree) ;
/*
   -----------------------
   read in the Tree object
   -----------------------
*/
tree = Tree_new() ;
Tree_readFromBinaryFile(tree, fp) ;
/*
   ---------------------
   read in the IV object
   ---------------------
*/
mapIV = IV_new() ;
IV_readFromBinaryFile(mapIV, fp) ;
/*
   ----------------------------
   initialize the DSTree object
   ----------------------------
*/
DSTree_init2(dstree, tree, mapIV) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to write an DSTree object to a file

   input --

      fn -- filename
        *.dstreeb -- binary
        *.dstreef -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ----------------------------------------------
*/
int
DSTree_writeToFile ( 
   DSTree   *dstree, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_writeToFile(%p,%s)"
    "\n bad input\n", dstree, fn) ; 
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
         fprintf(stderr, "\n error in DSTree_writeToFile(%p,%s)"
                 "\n unable to open file %s", dstree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DSTree_writeToBinaryFile(dstree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in DSTree_writeToFile(%p,%s)"
                 "\n unable to open file %s", dstree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DSTree_writeToFormattedFile(dstree, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in DSTree_writeToFile(%p,%s)"
                 "\n unable to open file %s", dstree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = DSTree_writeForHumanEye(dstree, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in DSTree_writeToFile(%p,%s)"
              "\n unable to open file %s", dstree, fn, fn) ;
      rc = 0 ;
   } else {
      rc = DSTree_writeForHumanEye(dstree, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to write an DSTree object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   --------------------------------------------------------
*/
int
DSTree_writeToFormattedFile ( 
   DSTree   *dstree, 
   FILE    *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fp == NULL || dstree->tree == NULL ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_writeToFormattedFile(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   exit(-1) ;
}
/*
   ---------------------------------
   write the Tree object to the file 
   ---------------------------------
*/
rc = Tree_writeToFormattedFile(dstree->tree, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing Tree to file\n", 
           dstree, fp, rc) ;
   return(0) ;
}
/*
   -------------------------------
   write the IV object to the file 
   -------------------------------
*/
rc = IV_writeToFormattedFile(dstree->mapIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in DSTree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing IV to file\n", 
           dstree, fp, rc) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to write an DSTree object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   -----------------------------------------------------
*/
int
DSTree_writeToBinaryFile ( 
   DSTree    *dstree, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fp == NULL || dstree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_writeToBinaryFile(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   exit(-1) ;
}
/*
   ---------------------------------
   write the Tree object to the file 
   ---------------------------------
*/
rc = Tree_writeToBinaryFile(dstree->tree, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in DSTree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing Tree to file\n", 
           dstree, fp, rc) ;
   return(0) ;
}
/*
   -------------------------------
   write the IV object to the file 
   -------------------------------
*/
rc = IV_writeToBinaryFile(dstree->mapIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in DSTree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing IV to file\n", 
           dstree, fp, rc) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to write an DSTree object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ----------------------------------------------------
*/
int
DSTree_writeForHumanEye ( 
   DSTree   *dstree, 
   FILE     *fp 
) {
int   rc ;

if ( dstree == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in DSTree_writeForHumanEye(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   exit(-1) ;
}
if ( (rc = DSTree_writeStats(dstree, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in DSTree_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from DSTree_writeStats(%p,%p)\n",
           dstree, fp, rc, dstree, fp) ;
   return(0) ;
}
fprintf(fp, "\n\n domain/separator tree") ;
Tree_writeForHumanEye(dstree->tree, fp) ;
fprintf(fp, "\n\n map from vertices to domains and separators") ;
IV_writeForHumanEye(dstree->mapIV, fp) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to write out the statistics for the DSTree object

   return value -- 1 if success, 0 otherwise

   created -- 96mar10, cca
   ------------------------------------------------------------
*/
int
DSTree_writeStats ( 
   DSTree   *dstree, 
   FILE     *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( dstree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in DSTree_writeStats(%p,%p)"
           "\n bad input\n", dstree, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, "\n DSTree : dstree object") ;
if ( rc < 0 ) { 
   fprintf(stderr, "\n fatal error in DSTree_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", dstree, fp, rc) ;
   return(0) ;
}
if ( dstree->tree != NULL && dstree->mapIV != NULL ) {
   rc = fprintf(fp, 
       "\n   %d domains and separators, %d vertices, occupies %d bytes",
       dstree->tree->n, IV_size(dstree->mapIV), DSTree_sizeOf(dstree)) ;
   if ( rc < 0 ) { 
      fprintf(stderr, "\n fatal error in DSTree_writeStats(%p,%p)"
              "\n rc = %d, return from fprintf\n", dstree, fp, rc) ;
      return(0) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
