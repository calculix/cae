/*  IO.c  */

#include "../BPG.h"

static const char *suffixb = ".bpgb" ;
static const char *suffixf = ".bpgf" ;

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to read in a BPG object from a file

   input --

      fn -- filename, must be *.bpgb or *.bpgf

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   ----------------------------------------------
*/
int
BPG_readFromFile ( 
   BPG    *bpg, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in BPG_readFromFile(%p,%s)"
           "\n bad input\n", bpg, fn) ;
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
         fprintf(stderr, "\n error in BPG_readFromFile(%p,%s)"
                 "\n unable to open file %s", bpg, fn, fn) ;
         rc = 0 ;
      } else {
         rc = BPG_readFromBinaryFile(bpg, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in BPG_readFromFile(%p,%s)"
                 "\n unable to open file %s", bpg, fn, fn) ;
         rc = 0 ;
      } else {
         rc = BPG_readFromFormattedFile(bpg, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in BPG_readFromFile(%p,%s)"
              "\n bad BPG file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              bpg, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in BPG_readFromFile(%p,%s)"
       "\n bad BPG file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       bpg, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read a BPG object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95oct06, cca
   --------------------------------------------------------
*/
int
BPG_readFromFormattedFile ( 
   BPG    *bpg, 
   FILE   *fp 
) {
Graph   *graph ;
int     nX, nY, rc ;
int     itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in BPG_readFromFormattedFile(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
BPG_clearData(bpg) ;
/*
   -----------------
   read in nX and nY
   -----------------
*/
if ( (rc = IVfscanf(fp, 2, itemp)) != 2 ) {
   fprintf(stderr, "\n error in BPG_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", bpg, fp, rc, 2) ;
   return(0) ;
}
nX = itemp[0] ;
nY = itemp[1] ;
/*
   --------------------------------------------------
   create the Graph IVL object, then read in its data
   --------------------------------------------------
*/
graph = Graph_new() ;
Graph_setDefaultFields(graph) ;
rc = Graph_readFromFormattedFile(graph, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in BPG_readFromFormattedFile(%p,%p)"
           "\n trying to read in Graph"
           "\n return code %d from Graph_readFormattedFile(%p,%p)",
           bpg, fp, rc, graph, fp) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
BPG_init(bpg, nX, nY, graph) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to read a BPG object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95oct06, cca
   ----------------------------------------------------
*/
int
BPG_readFromBinaryFile ( 
   BPG    *bpg, 
   FILE   *fp 
) {
Graph   *graph ;
int     nX, nY, rc ;
int     itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_readFromBinaryFile(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
BPG_clearData(bpg) ;
/*
   --------------------------------------------
   read in the two scalar parameters, nX and nY
   --------------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 2, fp)) != 2 ) {
   fprintf(stderr, "\n error in BPG_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", bpg, fp, rc, 2) ;
   return(0) ;
}
nX = itemp[0] ;
nY = itemp[1] ;
/*
   ----------------------------------------------
   create the Grapg object, then read in its data
   ----------------------------------------------
*/
graph = Graph_new() ;
Graph_setDefaultFields(graph) ;
rc = Graph_readFromBinaryFile(graph, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in BPG_readFromBinaryFile(%p,%p)"
           "\n trying to read in Graph"
           "\n return code %d from Graph_readBinaryFile(%p,%p)",
           bpg, fp, rc, graph, fp) ;
   return(0) ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
BPG_init(bpg, nX, nY, graph) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to write a BPG object to a file

   input --

      fn -- filename
        *.bpgb -- binary
        *.bpgf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   --------------------------------------------
*/
int
BPG_writeToFile ( 
   BPG    *bpg, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_writeToFile(%p,%s)"
    "\n bad input\n", bpg, fn) ; 
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
         fprintf(stderr, "\n error in BPG_writeToFile(%p,%s)"
                 "\n unable to open file %s", bpg, fn, fn) ;
         rc = 0 ;
      } else {
         rc = BPG_writeToBinaryFile(bpg, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in BPG_writeToFile(%p,%s)"
                 "\n unable to open file %s", bpg, fn, fn) ;
         rc = 0 ;
      } else {
         rc = BPG_writeToFormattedFile(bpg, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in BPG_writeToFile(%p,%s)"
                 "\n unable to open file %s", bpg, fn, fn) ;
         rc = 0 ;
      } else {
         rc = BPG_writeForHumanEye(bpg, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in BPG_writeToFile(%p,%s)"
              "\n unable to open file %s", bpg, fn, fn) ;
      rc = 0 ;
   } else {
      rc = BPG_writeForHumanEye(bpg, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write a BPG object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ------------------------------------------------------
*/
int
BPG_writeToFormattedFile ( 
   BPG    *bpg, 
   FILE   *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_writeToFormattedFile(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   return(0) ;
}
/*
   -----------------------------------
   write out the two scalar parameters
   -----------------------------------
*/
rc = fprintf(fp, "\n %d %d", bpg->nX, bpg->nY) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in BPG_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", bpg, fp, rc) ;
   return(0) ;
}
/*
   -------------------
   write out the graph
   -------------------
*/
rc = Graph_writeToFormattedFile(bpg->graph, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in BPG_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from Graph_writeToFormattedFile(%p,%p)" 
           "\n while attempting to write out graph\n",
           bpg, fp, rc, bpg->graph, fp) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write a BPG object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   ---------------------------------------------------
*/
int
BPG_writeToBinaryFile ( 
   BPG    *bpg, 
   FILE   *fp 
) {
int   rc ;
int   itemp[6] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_writeToBinaryFile(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   return(0) ;
}
/*
   -----------------------------------
   write out the two scalar parameters
   -----------------------------------
*/
itemp[0] = bpg->nX ;
itemp[1] = bpg->nY ;
rc = fwrite((void *) itemp, sizeof(int), 6, fp) ;
if ( rc != 6 ) {
   fprintf(stderr, "\n error in BPG_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", bpg, fp, rc, 6) ;
   return(0) ;
}
/*
   -------------------
   write out the graph
   -------------------
*/
rc = Graph_writeToBinaryFile(bpg->graph, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in BPG_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from Graph_writeToBinaryFile(%p,%p)" 
           "\n while attempting to write out graph\n",
           bpg, fp, rc, bpg->graph, fp) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write a BPG object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -------------------------------------------------
*/
int
BPG_writeForHumanEye ( 
   BPG    *bpg, 
   FILE   *fp 
) {
int   ierr, rc ;

if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in BPG_writeForHumanEye(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   exit(-1) ;
}
/*
   ------------------------
   write out the statistics
   ------------------------
*/
if ( (rc = BPG_writeStats(bpg, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in BPG_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from BPG_writeStats(%p,%p)\n",
           bpg, fp, rc, bpg, fp) ;
   return(0) ;
}
fflush(fp) ;
if ( bpg->graph != NULL ) {
/*
   --------------------------
   write out the Graph object
   --------------------------
*/
   fprintf(fp, "\n\n Graph object") ;
   rc = Graph_writeForHumanEye(bpg->graph, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, "\n fatal error in BPG_writeForHumanEye(%p,%p)"
              "\n rc = %d, return from Graph_writeForHumanEye(%p,%p)" 
              "\n while attempting to write out graph\n",
              bpg, fp, rc, bpg->graph, fp) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the BPG object

   return value -- 1 if success, 0 otherwise

   created -- 95oct06, cca
   -----------------------------------------------------------
*/
int
BPG_writeStats ( 
   BPG    *bpg, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in BPG_writeStats(%p,%p)"
           "\n bad input\n", bpg, fp) ;
   exit(-1) ;
}
if ( bpg->graph == NULL ) {
   fprintf(stderr, "\n warning in BPG_writeStats(%p,%p)"
           "\n bpg->graph = NULL\n", bpg, fp) ;
   return(1) ;
}
switch ( bpg->graph->type ) {
case 0 : 
   rc = fprintf(fp, "\n BPG : unweighted bpg object :") ;
   break ;
case 1 : 
   rc = fprintf(fp, "\n BPG : vertices weighted bpg object :") ;
   break ;
case 2 : 
   rc = fprintf(fp, "\n BPG : edges weighted bpg object :") ;
   break ;
case 3 : 
   rc = fprintf(fp, 
            "\n BPG : vertices and edges weighted bpg object :") ;
   break ;
default :
   fprintf(stderr, "\n fatal error in BPG_writeStats(%p,%p)"
           "\n invalid bpg->g->type = %d\n", 
           bpg, fp, bpg->graph->type) ;
   return(0) ;
}
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;
rc = fprintf(fp, " nX = %d, nY = %d", bpg->nX, bpg->nY) ;
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;
if ( bpg->graph != NULL ) {
   if ( bpg->graph->vwghts != NULL ) {
      rc = fprintf(fp, ", |X| = %d, |Y| = %d", 
                  IVsum(bpg->nX, bpg->graph->vwghts),
                  IVsum(bpg->nY, bpg->graph->vwghts + bpg->nX)) ;
   } else {
      rc = fprintf(fp, ", |X| = %d, |Y| = %d", bpg->nX, bpg->nY) ;
   }
}
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;

return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in BPG_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", bpg, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
