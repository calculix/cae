/*  IO.c  */

#include "../Graph.h"

#define BUFLEN 100000

static const char *suffixb = ".graphb" ;
static const char *suffixf = ".graphf" ;

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to read in a Graph object from a file

   input --

      fn -- filename, must be *.graphb or *.graphf

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   -------------------------------------------------
*/
int
Graph_readFromFile ( 
   Graph   *graph, 
   char    *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
           "\n bad input\n", graph, fn) ;
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
         fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
                 "\n unable to open file %s", graph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Graph_readFromBinaryFile(graph, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
                 "\n unable to open file %s", graph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Graph_readFromFormattedFile(graph, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
              "\n bad Graph file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              graph, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
       "\n bad Graph file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       graph, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to read in a Graph object from a CHACO file
 
   input --
 
      fn -- filename
 
   return value -- 1 if success, 0 if failure
 
   created -- 98sep20, jjs
   --------------------------------------------------------
*/
int
Graph_readFromChacoFile (
   Graph   *graph,
   char    *fn
) {
char    *rc ;
FILE    *fp;
int     nvtx, nedges, format;
char    string[BUFLEN], *s1, *s2;
int     k, v, vsize, w, vwghts, ewghts;
int     *adjncy, *weights, *vwghtsINT;
IVL     *adjIVL, *ewghtIVL;
/*
   ---------------
   check the input
   ---------------
*/
if ((graph == NULL) || (fn == NULL)) {
   fprintf(stderr, "\n error in Graph_readFromFile(%p,%s)"
           "\n bad input\n", graph, fn);
   return(0);
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Graph_clearData(graph);
/*
   ----------------------------------------------
   open file and read in nvtx, nedges, and format
   ----------------------------------------------
*/
if ((fp = fopen(fn, "r")) == (FILE*)NULL) {
   fprintf(stderr, "\n error in Graph_readFromChacoFile(%p,%s)"
           "\n unable to open file %s", graph, fn, fn);
   return(0);
}
/*
   -------------
   skip comments
   -------------
*/
do {
   rc = fgets(string, BUFLEN, fp) ;
   if ( rc == NULL ) {
      fprintf(stderr, "\n error in Graph_readFromChacoFile()"
             "\n error skipping comments in file %s\n", fn) ;
      return(0) ;
   }
} while ( string[0] == '%');
/*
   -------------------------------------------------
   read in # vertices, # edges and (optional) format
   -------------------------------------------------
*/
format = 0;
if (sscanf(string, "%d %d %d", &nvtx, &nedges, &format) < 2) {
   fprintf(stderr, "\n error in Graph_readFromChacoFile(%p,%s)"
           "\n unable to read header of file %s", graph, fn, fn);
   return(0);
}
ewghts = ((format % 10) > 0);
vwghts = (((format / 10) % 10) > 0);
if (format >= 100) {
   fprintf(stderr, "\n error in Graph_readFromChacoFile(%p,%s)"
           "\n unknown format", graph, fn);
   return(0);
}
/*
   ------------------------------------------------------------------
   initialize vector(s) to hold adjacency and (optional) edge weights
   ------------------------------------------------------------------
*/
adjncy = IVinit(nvtx, -1) ;
if ( ewghts ) {
   weights = IVinit(nvtx, -1) ;
} else {
   weights = NULL ;
}
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
nedges *= 2;
nedges += nvtx;
Graph_init1(graph, 2*ewghts+vwghts, nvtx, 0, nedges, 
            IVL_CHUNKED, IVL_CHUNKED);
adjIVL = graph->adjIVL;
if (ewghts) {
   ewghtIVL = graph->ewghtIVL;
   weights[0] = 0;                 /* self loops have no weight */
}
if (vwghts) vwghtsINT = graph->vwghts;
/*
   ---------------------------
   read in all adjacency lists
   ---------------------------
*/
k = 0;
for (v = 0; v < nvtx; v++) {
/*
   -------------
   skip comments
   -------------
*/
   do {
      rc = fgets(string, BUFLEN, fp);
      if ( rc == NULL ) {
         fprintf(stderr, "\n error in Graph_readFromChacoFile()"
                "\n error reading adjacency for vertex %d in file %s\n",
                v, fn) ;
         IVfree(adjncy) ;
         if ( weights != NULL ) {
            IVfree(weights) ;
         }
         return(0) ;
      }
   } while ( string[0] == '%');
/*
   -------------------------
   check for buffer overflow
   -------------------------
*/
   if (strlen(string) == BUFLEN-1) {
      fprintf(stderr, "\n error in Graph_readFromChacoFile(%p,%s)"
              "\n unable to read adjacency lists of file %s (line "
              "buffer too small)\n", graph, fn, fn);
      IVfree(adjncy) ;
      if ( weights != NULL ) {
         IVfree(weights) ;
      }
      return(0);
   }
/*
   ----------------------------------------------
   read in (optional) vertex weight, 
   adjacent vertices, and (optional) edge weights
   ----------------------------------------------
*/ 
   s1 = string;
   if (vwghts) vwghtsINT[v] = (int)strtol(string, &s1, 10);
   adjncy[0] = v;              /* insert self loop needed by spooles */
   if ( ewghts ) {
      weights[0] = 0;
   }
   vsize = 1;
   while ((w = (int)strtol(s1, &s2, 10)) > 0) {
      adjncy[vsize] = --w;     /* node numbering starts with 0 */
      s1 = s2;
      if (ewghts)
       { weights[vsize] = (int)strtol(s1, &s2, 10);
         s1 = s2;
       }
      vsize++;
   }
/*
   ---------------------------------
   sort the lists in ascending order
   ---------------------------------
*/
   if ( ewghts ) {
      IV2qsortUp(vsize, adjncy, weights) ;
   } else {
      IVqsortUp(vsize, adjncy) ;
   }
/*
   --------------------------------
   set the lists in the IVL objects
   --------------------------------
*/
   IVL_setList(adjIVL, v, vsize, adjncy);
   if (ewghts) IVL_setList(ewghtIVL, v, vsize, weights);
   k += vsize;
}
/*
   -----------------------------------
   close the file and do a final check
   -----------------------------------
*/
fclose(fp);
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(adjncy) ;
if ( weights != NULL ) {
   IVfree(weights) ;
}
/*
   ----------------
   check for errors
   ----------------
*/
if ((k != nedges) || (v != nvtx)) {
   fprintf(stderr, "\n error in Graph_readFromChacoFile()"
           "\n number of nodes/edges does not match with header of %s"
           "\n k %d, nedges %d, v %d, nvtx %d\n", 
           fn, k, nedges, v, nvtx);
   return(0);
}
return(1); }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read a Graph object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95sep29, cca
   --------------------------------------------------------
*/
int
Graph_readFromFormattedFile ( 
   Graph   *graph, 
   FILE    *fp 
) {
int   nedges, nvbnd, nvtx, rc, totewght, totvwght, type ;
int   itemp[6] ;
int   *vwghts ;
IVL   *adjIVL, *ewghtIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
           "\n bad input\n", graph, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Graph_clearData(graph) ;
/*
   ---------------------------------------------
   read in the six scalar parameters
   type, nvtx, nvbnd, nedges, totvwght, totewght
   ---------------------------------------------
*/
if ( (rc = IVfscanf(fp, 6, itemp)) != 6 ) {
   fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", graph, fp, rc, 6) ;
   return(0) ;
}
type     = itemp[0] ;
nvtx     = itemp[1] ;
nvbnd    = itemp[2] ;
nedges   = itemp[3] ;
totvwght = itemp[4] ;
totewght = itemp[5] ;
/*
   --------------------------------------------------
   create the adjIVL IVL object, 
   set its type to IVL_CHUNKED, then read in its data
   --------------------------------------------------
*/
adjIVL = IVL_new() ;
IVL_setDefaultFields(adjIVL) ;
adjIVL->type = IVL_CHUNKED ;
rc = IVL_readFromFormattedFile(adjIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
           "\n trying to read in adjIVL"
           "\n return code %d from IVL_readFormattedFile(%p,%p)",
           graph, fp, rc, adjIVL, fp) ;
   return(0) ;
}
if ( type % 2 == 1 ) {
   int   nvtot, wght ;
/*
   --------------------------
   vertex weights are present
   --------------------------
*/
   nvtot = nvtx + nvbnd ;
   vwghts = IVinit2(nvtot) ;
   if ( (rc = IVfscanf(fp, nvtot, vwghts)) != nvtot ) {
      fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
              "\n %d items of %d read\n", graph, fp, rc, nvtot) ;
      return(0) ;
   }
   wght = IVsum(nvtot, vwghts) ;
   if ( wght != totvwght ) {
      fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
              "\n totvwght = %d, IVsum(vwghts) = %d\n",
              graph, fp, totvwght, wght) ;
      return(0) ;
   }
} else {
   vwghts = NULL ;
}
if ( type >= 2 ) {
   int   wght ;
/*
   -----------------------------------------------------
   edge weights are present, create the ewghtIVL object, 
   set its type to IVL_CHUNKED, then read in its data
   -----------------------------------------------------
*/
   ewghtIVL = IVL_new() ;
   IVL_setDefaultFields(ewghtIVL) ;
   ewghtIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromFormattedFile(ewghtIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
              "\n trying to read in ewghtIVL"
              "\n return code %d from IVL_readFormattedFile(%p,%p)",
              graph, fp, rc, ewghtIVL, fp) ;
      return(0) ;
   }
   wght = IVL_sum(ewghtIVL) ;
   if ( wght != totewght ) {
      fprintf(stderr, "\n error in Graph_readFromFormattedFile(%p,%p)"
              "\n totewght = %d, IVL_sum(ewghtIVL) = %d\n",
              graph, fp, totewght, wght) ;
      return(0) ;
   }
} else {
   ewghtIVL = NULL ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Graph_init2(graph, type, nvtx, nvbnd, nedges, totvwght, totewght,
            adjIVL, vwghts, ewghtIVL) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to read a Graph object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95sep29, cca
   ----------------------------------------------------
*/
int
Graph_readFromBinaryFile ( 
   Graph   *graph, 
   FILE    *fp 
) {
int   nedges, nvbnd, nvtx, rc, totewght, totvwght, type ;
int   itemp[6] ;
int   *vwghts ;
IVL   *adjIVL, *ewghtIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_readFromBinaryFile(%p,%p)"
           "\n bad input\n", graph, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
Graph_clearData(graph) ;
/*
   ---------------------------------------------
   read in the six scalar parameters
   type, nvtx, nvbnd, nedges, totvwght, totewght
   ---------------------------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 6, fp)) != 6 ) {
   fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", graph, fp, rc, 6) ;
   return(0) ;
}
type     = itemp[0] ;
nvtx     = itemp[1] ;
nvbnd    = itemp[2] ;
nedges   = itemp[3] ;
totvwght = itemp[4] ;
totewght = itemp[5] ;
/*
   --------------------------------------------------
   create the adjIVL IVL object, 
   set its type to IVL_CHUNKED, then read in its data
   --------------------------------------------------
*/
adjIVL = IVL_new() ;
IVL_setDefaultFields(adjIVL) ;
adjIVL->type = IVL_CHUNKED ;
rc = IVL_readFromBinaryFile(adjIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
           "\n trying to read in adjIVL"
           "\n return code %d from IVL_readBinaryFile(%p,%p)",
           graph, fp, rc, adjIVL, fp) ;
   return(0) ;
}
if ( type % 2 == 1 ) {
   int   nvtot, wght ;
/*
   --------------------------
   vertex weights are present
   --------------------------
*/
   nvtot  = nvtx + nvbnd ;
   vwghts = IVinit2(nvtot) ;
   if ( (rc = fread((void *) vwghts, sizeof(int), nvtot, fp)) != nvtot){
      fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", graph, fp, rc, nvtx+nvbnd) ;
      return(0) ;
   }
   wght = IVsum(nvtot, vwghts) ;
   if ( wght != totvwght ) {
      fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
              "\n totvwght = %d, IVsum(vwghts) = %d\n",
              graph, fp, totvwght, wght) ;
      return(0) ;
   }
} else {
   vwghts = NULL ;
}
if ( type > 2 ) {
   int   wght ;
/*
   -----------------------------------------------------
   edge weights are present, create the ewghtIVL object, 
   set its type to IVL_CHUNKED, then read in its data
   -----------------------------------------------------
*/
   ewghtIVL = IVL_new() ;
   IVL_setDefaultFields(ewghtIVL) ;
   ewghtIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromBinaryFile(ewghtIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
              "\n trying to read in ewghtIVL"
              "\n return code %d from IVL_readBinaryFile(%p,%p)",
              graph, fp, rc, ewghtIVL, fp) ;
      return(0) ;
   }
   wght = IVL_sum(ewghtIVL) ;
   if ( wght != totewght ) {
      fprintf(stderr, "\n error in Graph_readFromBinaryFile(%p,%p)"
              "\n totewght = %d, IVL_sum(ewghtIVL) = %d\n",
              graph, fp, totewght, wght) ;
      return(0) ;
   }
} else {
   ewghtIVL = NULL ;
}
/*
   ---------------------
   initialize the object
   ---------------------
*/
Graph_init2(graph, type, nvtx, nvbnd, nedges, totvwght, totewght,
            adjIVL, vwghts, ewghtIVL) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to write a Graph object to a file

   input --

      fn -- filename
        *.graphb -- binary
        *.graphf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   --------------------------------------------
*/
int
Graph_writeToFile ( 
   Graph   *graph, 
   char    *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFile(%p,%s)"
    "\n bad input\n", graph, fn) ; 
   return(0) ;
}
if ( graph->type < 0 || 3 < graph->type ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFile(%p,%s)"
           "\n bad type = %d", graph, fn, graph->type) ;
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
         fprintf(stderr, "\n error in Graph_writeToFile(%p,%s)"
                 "\n unable to open file %s", graph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Graph_writeToBinaryFile(graph, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in Graph_writeToFile(%p,%s)"
                 "\n unable to open file %s", graph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Graph_writeToFormattedFile(graph, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in Graph_writeToFile(%p,%s)"
                 "\n unable to open file %s", graph, fn, fn) ;
         rc = 0 ;
      } else {
         rc = Graph_writeForHumanEye(graph, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in Graph_writeToFile(%p,%s)"
              "\n unable to open file %s", graph, fn, fn) ;
      rc = 0 ;
   } else {
      rc = Graph_writeForHumanEye(graph, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write a Graph object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   ------------------------------------------------------
*/
int
Graph_writeToFormattedFile ( 
   Graph   *graph, 
   FILE    *fp 
) {
int   ierr, rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
           "\n bad input\n", graph, fp) ;
   return(0) ;
}
if ( graph->type < 0 || 3 < graph->type ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
           "\n bad type = %d", graph, fp, graph->type) ;
   return(0) ;
}
/*
   -----------------------------------
   write out the six scalar parameters
   -----------------------------------
*/
rc = fprintf(fp, "\n %d %d %d %d %d %d", 
             graph->type, graph->nvtx, graph->nvbnd,
             graph->nedges, graph->totvwght, graph->totewght) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", graph, fp, rc) ;
   return(0) ;
}
/*
   ---------------------------------
   write out the adjacency structure
   ---------------------------------
*/
rc = IVL_writeToFormattedFile(graph->adjIVL, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from IVL_writeToFormattedFile(%p,%p)" 
           "\n while attempting to write out adjIVL\n",
           graph, fp, rc, graph->adjIVL, fp) ;
   return(0) ;
}
/*
   -----------------------------------------
   write out the vwghts[] vector, if present
   -----------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   if ( graph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
              "\n graph->type = %d, graph->vwghts == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   IVfp80(fp, graph->nvtx+graph->nvbnd, graph->vwghts, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
              "\n ierr = %d, return from vwghts[] IVfp80\n", 
              graph, fp, ierr) ;
      return(0) ;
   }
}
/*
   -------------------------------------
   write out the edge weights IVL object 
   -------------------------------------
*/
if ( graph->type >= 2 ) {
   if ( graph->ewghtIVL == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
              "\n graph->type = %d, graph->ewghtIVL == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   rc = IVL_writeToFormattedFile(graph->ewghtIVL, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToFormattedFile(%p,%p)"
              "\n rc = %d, return from IVL_writeToFormattedFile(%p,%p)"
              "\n while attempting to write out ewghtIVL\n",
              graph, fp, rc, graph->ewghtIVL, fp) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write a Graph object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   ---------------------------------------------------
*/
int
Graph_writeToBinaryFile ( 
   Graph    *graph, 
   FILE   *fp 
) {
int   rc ;
int   itemp[6] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
           "\n bad input\n", graph, fp) ;
   return(0) ;
}
if ( graph->type < 0 || 3 < graph->type ) {
   fprintf(stderr, "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
           "\n bad type = %d", graph, fp, graph->type) ;
   return(0) ;
}
/*
   -----------------------------------
   write out the six scalar parameters
   -----------------------------------
*/
itemp[0] = graph->type     ;
itemp[1] = graph->nvtx     ;
itemp[2] = graph->nvbnd    ;
itemp[3] = graph->nedges   ;
itemp[4] = graph->totvwght ;
itemp[5] = graph->totewght ;
rc = fwrite((void *) itemp, sizeof(int), 6, fp) ;
if ( rc != 6 ) {
   fprintf(stderr, "\n error in Graph_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", graph, fp, rc, 6) ;
   return(0) ;
}
/*
   ---------------------------------
   write out the adjacency structure
   ---------------------------------
*/
rc = IVL_writeToBinaryFile(graph->adjIVL, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from IVL_writeToBinaryFile(%p,%p)" 
           "\n while attempting to write out adjIVL\n",
           graph, fp, rc, graph->adjIVL, fp) ;
   return(0) ;
}
/*
   -----------------------------------------
   write out the vwghts[] vector, if present
   -----------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   if ( graph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
              "\n graph->type = %d, graph->vwghts == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   rc = fwrite((void *) graph->vwghts, sizeof(int), 
               graph->nvtx + graph->nvbnd, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from vwghts[] fwrite\n", 
              graph, fp, rc) ;
      return(0) ;
   }
}
/*
   -------------------------------------
   write out the edge weights IVL object 
   -------------------------------------
*/
if ( graph->type >= 2 ) {
   if ( graph->ewghtIVL == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
              "\n graph->type = %d, graph->ewghtIVL == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   rc = IVL_writeToBinaryFile(graph->ewghtIVL, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeToBinaryFile(%p,%p)"
              "\n rc = %d, return from IVL_writeToBinaryFile(%p,%p)"
              "\n while attempting to write out ewghtIVL\n",
              graph, fp, rc, graph->ewghtIVL, fp) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write a Graph object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -------------------------------------------------
*/
int
Graph_writeForHumanEye ( 
   Graph    *graph, 
   FILE   *fp 
) {
int   ierr, rc ;

if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_writeForHumanEye(%p,%p)"
           "\n bad input\n", graph, fp) ;
   exit(-1) ;
}
/*
   ------------------------
   write out the statistics
   ------------------------
*/
if ( (rc = Graph_writeStats(graph, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in Graph_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from Graph_writeStats(%p,%p)\n",
           graph, fp, rc, graph, fp) ;
   return(0) ;
}
if ( graph->adjIVL != NULL ) {
/*
   ----------------------------------
   write out the adjacency IVL object
   ----------------------------------
*/
   fprintf(fp, "\n\n adjacency IVL object") ;
   rc = IVL_writeForHumanEye(graph->adjIVL, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, "\n fatal error in Graph_writeForHumanEye(%p,%p)"
              "\n rc = %d, return from IVL_writeForHumanEye(%p,%p)" 
              "\n while attempting to write out adjIVL\n",
              graph, fp, rc, graph->adjIVL, fp) ;
      return(0) ;
   }
}
/*
   ----------------------------------------------
   write out the vertex weights vector if present
   ----------------------------------------------
*/
if ( graph->type % 2 == 1 ) {
   if ( graph->vwghts == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeForHumanEye(%p,%p)"
              "\n graph->type = %d, graph->vwghts == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   fprintf(fp, "\n\n vertex weights ") ;
   IVfp80(fp, graph->nvtx + graph->nvbnd, graph->vwghts, 80, &ierr) ;
   if ( ierr < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeForHumanEye(%p,%p)"
              "\n ierr = %d, return from vwghts[] IVfp80\n", 
              graph, fp, ierr) ;
      return(0) ;
   }
}
/*
   ------------------------------------------------
   write out the edge weight IVL object, if present
   ------------------------------------------------
*/
if ( graph->type >= 2 ) {
   if ( graph->ewghtIVL == NULL ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeForHumanEye(%p,%p)"
              "\n graph->type = %d, graph->ewghtIVL == NULL\n",
              graph, fp, graph->type) ;
      return(0) ;
   }
   fprintf(fp, "\n\n edge weights IVL object") ;
   rc = IVL_writeForHumanEye(graph->ewghtIVL, fp) ;
   if ( rc < 0 ) {
      fprintf(stderr, 
              "\n fatal error in Graph_writeForHumanEye(%p,%p)"
              "\n rc = %d, return from IVL_writeForHumanEye(%p,%p)"
              "\n while attempting to write out ewghtIVL\n",
              graph, fp, rc, graph->ewghtIVL, fp) ;
      return(0) ;
   }
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the Graph object

   return value -- 1 if success, 0 otherwise

   created -- 95sep29, cca
   -----------------------------------------------------------
*/
int
Graph_writeStats ( 
   Graph    *graph, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( graph == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in Graph_writeStats(%p,%p)"
           "\n bad input\n", graph, fp) ;
   exit(-1) ;
}
switch ( graph->type ) {
case 0 : 
   rc = fprintf(fp, "\n Graph : unweighted graph object :") ;
   break ;
case 1 : 
   rc = fprintf(fp, "\n Graph : vertices weighted graph object :") ;
   break ;
case 2 : 
   rc = fprintf(fp, "\n Graph : edges weighted graph object :") ;
   break ;
case 3 : 
   rc = fprintf(fp, 
            "\n Graph : vertices and edges weighted graph object :") ;
   break ;
default :
   fprintf(stderr, "\n fatal error in Graph_writeStats(%p,%p)"
           "\n invalid graph->type = %d\n", graph, fp, graph->type) ;
   return(0) ;
}
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;
rc = fprintf(fp, 
             "\n %d internal vertices, %d boundary vertices, %d edges",
             graph->nvtx, graph->nvbnd, graph->nedges) ;
if ( rc < 0 ) { goto IO_error ; }
fflush(fp) ;
rc = fprintf(fp, 
             "\n %d internal vertex weight, %d boundary vertex weight",
(graph->vwghts != NULL) ? IVsum(graph->nvtx, graph->vwghts) : graph->nvtx,
(graph->vwghts != NULL) ? IVsum(graph->nvbnd, graph->vwghts + graph->nvtx) : graph->nvbnd) ;
if ( rc < 0 ) { goto IO_error ; }
if ( graph->type >= 2 ) {
   rc = fprintf(fp, "\n %d total edge weight",
                graph->totewght) ;
}
if ( rc < 0 ) { goto IO_error ; }

return(1) ;

IO_error :
   fprintf(stderr, "\n fatal error in Graph_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", graph, fp, rc) ;
   return(0) ;
}

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   write out a graph to a METIS file

   created -- 95oct18, cca
   ---------------------------------
*/
int
Graph_writeToMetisFile (
   Graph   *g,
   FILE    *fp
) {
int   ii, nedge, nvtx, v, vsize, w ;
int   *vewghts, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Graph_writeToMetisFile(%p,%p)"
           "\n bad input\n", g, fp) ;
   exit(-1) ;
}
nvtx = g->nvtx ;
nedge = (g->nedges - nvtx)/2 ;
switch ( g->type ) {
case 0 : 
   fprintf(fp, " %d %d   ", nvtx, nedge) ; 
   for ( v = 0 ; v < nvtx ; v++ ) {
      fprintf(fp, "\n ") ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w != v && w < nvtx ) {
            fprintf(fp, " %d", w + 1) ;
         }
      }
   }
   break ;
case 1 : 
   fprintf(fp, " %d %d 10", nvtx, nedge) ; 
   for ( v = 0 ; v < nvtx ; v++ ) {
      fprintf(fp, "\n %d", g->vwghts[v]) ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w != v && w < nvtx ) {
            fprintf(fp, " %d", w + 1) ;
         }
      }
   }
   break ;
case 2 : 
   fprintf(fp, " %d %d  1", nvtx, nedge) ; 
   for ( v = 0 ; v < nvtx ; v++ ) {
      fprintf(fp, "\n") ;
      Graph_adjAndEweights(g, v, &vsize, &vadj, &vewghts) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w != v && w < nvtx ) {
            fprintf(fp, " %d %d", w + 1, vewghts[ii]) ;
         }
      }
   }
   break ;
case 3 : 
   fprintf(fp, " %d %d 11", nvtx, nedge) ; 
   for ( v = 0 ; v < nvtx ; v++ ) {
      fprintf(fp, "\n %d", g->vwghts[v]) ;
      Graph_adjAndEweights(g, v, &vsize, &vadj, &vewghts) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w != v && w < nvtx ) {
            fprintf(fp, " %d %d", w + 1, vewghts[ii]) ;
         }
      }
   }
   break ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
