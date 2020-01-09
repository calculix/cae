/*  IO.c  */

#include "../Tree.h"

static const char *suffixb = ".treeb" ;
static const char *suffixf = ".treef" ;

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- to read in an Tree object from a file

   input --

      fn -- filename, must be *.treeb or *.treef

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   ------------------------------------------------
*/
int
Tree_readFromFile ( 
   Tree   *tree, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in Tree_readFromFile(%p,%s)"
           "\n bad input\n", tree, fn) ;
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
         fprintf(stderr, "\n error in Tree_readFromFile(%p,%s)"
                 "\n unable to open file %s", tree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Tree_readFromBinaryFile(tree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in Tree_readFromFile(%p,%s)"
                 "\n unable to open file %s", tree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Tree_readFromFormattedFile(tree, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in Tree_readFromFile(%p,%s)"
              "\n bad Tree file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              tree, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in Tree_readFromFile(%p,%s)"
       "\n bad Tree file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       tree, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to read an Tree object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   -------------------------------------------------------
*/
int
Tree_readFromFormattedFile ( 
   Tree   *tree, 
   FILE   *fp 
) {
int   rc ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Tree_readFromFormattedFile(%p,%p)"
           "\n bad input\n", tree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Tree_clearData(tree) ;
/*
   ---------------------------------
   read in the two scalar parameters
   number of nodes and the root
   ---------------------------------
*/
if ( (rc = IVfscanf(fp, 2, itemp)) != 2 ) {
   fprintf(stderr, "\n error in Tree_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", tree, fp, rc, 2) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Tree_init1(tree, itemp[0]) ;
tree->root = itemp[1] ;
/*
   -----------------------
   now read in the indices
   -----------------------
*/
if ( (rc = IVfscanf(fp, tree->n, tree->par)) != tree->n ) {
   fprintf(stderr, "\n par: error in Tree_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}
if ( (rc = IVfscanf(fp, tree->n, tree->fch)) != tree->n ) {
   fprintf(stderr, "\n fch: error in Tree_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}
if ( (rc = IVfscanf(fp, tree->n, tree->sib)) != tree->n ) {
   fprintf(stderr, "\n sib: error in Tree_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to read an Tree object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95nov15, cca
   ----------------------------------------------------
*/
int
Tree_readFromBinaryFile ( 
   Tree    *tree, 
   FILE   *fp 
) {
int   rc ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_readFromBinaryFile(%p,%p)"
           "\n bad input\n", tree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Tree_clearData(tree) ;
/*
   ---------------------------------
   read in the two scalar parameters
   number of nodes and the root
   ---------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 2, fp)) != 2 ) {
   fprintf(stderr, "\n error in Tree_readFromBinaryFile(%p,%p)"
           "\n itemp(2) : %d items of %d read\n", tree, fp, rc, 2) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Tree_init1(tree, itemp[0]) ;
tree->root = itemp[1] ;
/*
   -------------------
   read in the vectors
   -------------------
*/
if ( (rc = fread((void *) tree->par, sizeof(int), tree->n, fp)) 
     != tree->n ) {
   fprintf(stderr, "\n par : error in Tree_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}
if ( (rc = fread((void *) tree->fch, sizeof(int), tree->n, fp)) 
     != tree->n ) {
   fprintf(stderr, "\n fch : error in Tree_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}
if ( (rc = fread((void *) tree->sib, sizeof(int), tree->n, fp)) 
     != tree->n ) {
   fprintf(stderr, "\n sib : error in Tree_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", 
           tree, fp, rc, tree->n) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to write an Tree object to a file

   input --

      fn -- filename
        *.treeb -- binary
        *.treef -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   --------------------------------------------
*/
int
Tree_writeToFile ( 
   Tree   *tree, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_writeToFile(%p,%s)"
    "\n bad input\n", tree, fn) ; 
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
         fprintf(stderr, "\n error in Tree_writeToFile(%p,%s)"
                 "\n unable to open file %s", tree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Tree_writeToBinaryFile(tree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in Tree_writeToFile(%p,%s)"
                 "\n unable to open file %s", tree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Tree_writeToFormattedFile(tree, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in Tree_writeToFile(%p,%s)"
                 "\n unable to open file %s", tree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Tree_writeForHumanEye(tree, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in Tree_writeToFile(%p,%s)"
              "\n unable to open file %s", tree, fn, fn) ;
      rc = 0 ;
   } else {
      rc = Tree_writeForHumanEye(tree, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write an Tree object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ------------------------------------------------------
*/
int
Tree_writeToFormattedFile ( 
   Tree   *tree, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fp == NULL || tree->n <= 0 ) {
   fprintf(stderr, "\n fatal error in Tree_writeToFormattedFile(%p,%p)"
           "\n bad input\n", tree, fp) ;
   exit(-1) ;
}
/*
   -------------------------------------
   write out the two scalar parameters
   -------------------------------------
*/
rc = fprintf(fp, "\n %d %d", tree->n, tree->root) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in Tree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", tree, fp, rc) ;
   return(0) ;
}
/*
   ---------------------
   write out the vectors
   ---------------------
*/
IVfp80(fp, tree->n, tree->par, 80, &ierr) ;
IVfp80(fp, tree->n, tree->fch, 80, &ierr) ;
IVfp80(fp, tree->n, tree->sib, 80, &ierr) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write an Tree object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
Tree_writeToBinaryFile ( 
   Tree    *tree, 
   FILE   *fp 
) {
int   rc ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fp == NULL || tree->n <= 0 ) {
   fprintf(stderr, "\n fatal error in Tree_writeToBinaryFile(%p,%p)"
           "\n bad input\n", tree, fp) ;
   exit(-1) ;
}
itemp[0] = tree->n    ;
itemp[1] = tree->root ;
rc = fwrite((void *) itemp, sizeof(int), 2, fp) ;
if ( rc != 2 ) {
   fprintf(stderr, "\n error in Tree_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", tree, fp, rc, 2) ;
   return(0) ;
}
rc = fwrite((void *) tree->par, sizeof(int), tree->n, fp) ;
if ( rc != tree->n ) {
   fprintf(stderr, "\n error in Tree_writeToBinaryFile(%p,%p)"
           "\n tree->par, %d of %d items written\n",
           tree, fp, rc, tree->n) ;
   return(0) ;
}
rc = fwrite((void *) tree->fch, sizeof(int), tree->n, fp) ;
if ( rc != tree->n ) {
   fprintf(stderr, "\n error in Tree_writeToBinaryFile(%p,%p)"
           "\n tree->fch, %d of %d items written\n",
           tree, fp, rc, tree->n) ;
   return(0) ;
}
rc = fwrite((void *) tree->sib, sizeof(int), tree->n, fp) ;
if ( rc != tree->n ) {
   fprintf(stderr, "\n error in Tree_writeToBinaryFile(%p,%p)"
           "\n tree->sib, %d of %d items written\n",
           tree, fp, rc, tree->n) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- to write an Tree object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   --------------------------------------------------
*/
int
Tree_writeForHumanEye ( 
   Tree    *tree, 
   FILE   *fp 
) {
int   rc, v ;

if ( tree == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Tree_writeForHumanEye(%p,%p)"
           "\n bad input\n", tree, fp) ;
   exit(-1) ;
}
if ( (rc = Tree_writeStats(tree, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in Tree_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from Tree_writeStats(%p,%p)\n",
           tree, fp, rc, tree, fp) ;
   return(0) ;
}
fprintf(fp, "\n vertex   parent   fchild   sibling") ;
for ( v = 0 ; v < tree->n ; v++ ) {
   fprintf(fp, "\n %5d %9d %9d %9d :", 
           v, tree->par[v], tree->fch[v], tree->sib[v]) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   purpose -- to write out the statistics for the Tree object

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ----------------------------------------------------------
*/
int
Tree_writeStats ( 
   Tree    *tree, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Tree_writeStats(%p,%p)"
           "\n bad input\n", tree, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, 
        "\n Tree : tree object, %d vertices, root = %d, takes %d bytes",
        tree->n, tree->root, Tree_sizeOf(tree)) ;
if ( rc < 0 ) { 
   fprintf(stderr, "\n fatal error in Tree_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", tree, fp, rc) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
