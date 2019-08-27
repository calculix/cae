/*  makeGraph.c  */

#include "../InpMtx.h"
#include "../../Graph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------
   read in an adjacency structure file 
   and create a unit weight Graph object

   nvtx nadj
   offsets[nvtx+1]
   indices[nadj]

   created -- 96nov15, cca
   -------------------------------------
*/
{
InpMtx   *inpmtx ;
double   t1, t2 ;
int      flag, ierr, ii, msglvl, nadj, nvtx, rc, v, vsize ;
int      *adjncy, *offsets, *vadj ;
Graph    *graph ;
FILE     *fp, *msgFile ;

if ( argc != 6 ) {
   fprintf(stdout, 
    "\n\n usage : %s msglvl msgFile inAdjacencyFile outGraphFile flag"
    "\n    msglvl          -- message level"
    "\n    msgFile         -- message file"
    "\n    inAdjacencyFile -- input file"
    "\n    outGraphFile    -- output file, must be *.graphf or *.graphb"
    "\n    flag            -- flag for 0-based or 1-based addressing"
      "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
flag = atoi(argv[5]) ;
if ( flag < 0 || flag > 1 ) {
   fprintf(stderr, "\n fatal error in makeGraph %d %s %s %s %d"
   "\n flag must be 0 if indices and offsets are zero-based like C"
   "\n flag must be 1 if indices and offsets are one-based like Fortran"
   "\n", msglvl, argv[2], argv[3], argv[4], flag) ;
   exit(-1) ;
}
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n flag     -- %d" 
        "\n",
        argv[0], msglvl, argv[2], argv[3], argv[4], flag) ;
fflush(msgFile) ;
/*
   -------------------
   open the input file
   -------------------
*/
if ( (fp = fopen(argv[3], "r")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n", argv[0], argv[3]) ;
   return(-1) ;
}
/*
   -------------------------------------------------
   read in the number of vertices and adjacency size
   -------------------------------------------------
*/
if ( (rc = fscanf(fp, " %d %d", &nvtx, &nadj)) != 2 ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n %d of %d items read\n", argv[0], rc, 2) ;
   exit(-1) ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n adjacency file has %d vertices and %d indices",
           nvtx, nadj) ;
   fflush(msgFile) ;
}
/*
   -------------------
   read in the offsets
   -------------------
*/
offsets = IVinit(nvtx + 1, -1) ;
rc = IVfscanf(fp, nvtx + 1, offsets) ;
if ( rc != nvtx + 1 ) {
   fprintf(stderr, "\n fatal error in %s reading offsets"
           "\n %d of %d items read\n", 
           argv[0], rc, nvtx + 1) ;
   exit(-1) ;
}
if ( flag == 1 ) {
/*
   ---------------------------------------------
   decrement the offsets by one for C addressing
   ---------------------------------------------
*/
   for ( ii = 0 ; ii < nvtx + 1 ; ii++ ) {
      offsets[ii]-- ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n offsets") ;
   IVfp80(msgFile, nvtx+1, offsets, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   read in the adjacency
   ---------------------
*/
adjncy = IVinit(nadj, -1) ;
rc = IVfscanf(fp, nadj, adjncy) ;
if ( rc != nadj ) {
   fprintf(stderr, "\n fatal error in %s reading adjncy"
           "\n %d of %d items read\n", 
           argv[0], rc, nadj) ;
   exit(-1) ;
}
if ( flag == 1 ) {
/*
   ---------------------------------------------
   decrement the offsets by one for C addressing
   ---------------------------------------------
*/
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      adjncy[ii]-- ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n adjncy") ;
   IVfp80(msgFile, nadj, adjncy, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------
   check the adjacency structure for valid entries
   -----------------------------------------------
*/
for ( ii = 0 ; ii < nadj ; ii++ ) {
   v = adjncy[ii] ;
   if ( v < 0 || v >= nvtx ) {
      fprintf(stderr, "\n fatal error in mkGraph"
              "\n adjncy[%d] = %d, out of range in [0,%d]"
              "\n ", ii, v, nvtx-1) ;
      exit(-1) ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n adjncy") ;
   IVfp80(msgFile, nadj, adjncy, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   set up the InpMtx object
   -------------------------
*/
inpmtx = InpMtx_new() ;
InpMtx_init(inpmtx, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, 
            nvtx + 2*nadj, nvtx) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   InpMtx_inputEntry(inpmtx, v, v) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n inputting entry (%d,%d)", v, v) ;
      fflush(msgFile) ;
   }
   vsize = offsets[v+1] - offsets[v] ;
   vadj  = &adjncy[offsets[v]] ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n inputting row %d :", v) ;
      IVfp80(msgFile, vsize, vadj, 22, &ierr) ;
      fflush(msgFile) ;
   }
   InpMtx_inputRow(inpmtx, v, vsize, vadj) ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n inputting column %d :", v) ;
      IVfp80(msgFile, vsize, vadj, 22, &ierr) ;
      fflush(msgFile) ;
   }
   InpMtx_inputColumn(inpmtx, v, vsize, vadj) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n raw data") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n after changing to vector storage mode") ;
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
}
/*
   -----------------------------------
   free the offsets and adjncy vectors
   -----------------------------------
*/
IVfree(offsets) ;
IVfree(adjncy) ;
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
offsets = IVinit(nvtx + 1, -1) ;
IVcopy(nvtx, offsets, IV_entries(&inpmtx->offsetsIV)) ;
offsets[nvtx] = inpmtx->nent ;
adjncy = IV_entries(&inpmtx->ivec2IV) ;
graph  = Graph_new() ;
Graph_fillFromOffsets(graph, nvtx, offsets, adjncy, 0) ;
graph->type = 0 ;
graph->totvwght = nvtx ;
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n graph has %d vertices and %d edges",
           graph->nvtx, graph->nedges) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
}
/*
   --------------------------
   write out the Graph object
   --------------------------
*/
if ( strcmp(argv[4], "none") != 0 ) {
   MARKTIME(t1) ;
   rc = Graph_writeToFile(graph, argv[4]) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, argv[4]) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_writeToFile(%p,%s)",
           rc, graph, argv[4]) ;
}
/*
   ---------------------------------
   free the Graph and InpMtx object
   ---------------------------------
*/
Graph_free(graph) ;
InpMtx_free(inpmtx) ;
IVfree(offsets) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
