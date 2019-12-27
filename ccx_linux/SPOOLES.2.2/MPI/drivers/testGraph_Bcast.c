/*  testGraph_Bcast.c  */

#include "../spoolesMPI.h"
#include "../../Drand.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------------------------
   this program tests the Graph_MPI_Bcast() method

   (1) process root generates a random Graph object
       and computes its checksum
   (2) process root broadcasts the Graph object to the other processors
   (3) each process computes the checksum of its Graph object
   (4) the checksums are compared on root

   created -- 98sep10, cca
   --------------------------------------------------------------------
*/
{
char         *buffer ;
double       chksum, t1, t2 ;
double       *sums ;
Drand        drand ;
int          iproc, length, loc, msglvl, myid, nitem, nproc, 
             nvtx, root, seed, size, type, v ;
int          *list ;
FILE         *msgFile ;
Graph        *graph ;
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
fprintf(stdout, "\n process %d of %d, argc = %d", myid, nproc, argc) ;
fflush(stdout) ;
if ( argc != 8 ) {
   fprintf(stdout, 
           "\n\n usage : %s msglvl msgFile type nvtx nitem root seed "
           "\n    msglvl      -- message level"
           "\n    msgFile     -- message file"
           "\n    type        -- type of graph"
           "\n    nvtx        -- # of vertices"
           "\n    nitem       -- # of items used to generate graph"
           "\n    root        -- root processor for broadcast"
           "\n    seed        -- random number seed"
           "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else {
   length = strlen(argv[2]) + 1 + 4 ;
   buffer = CVinit(length, '\0') ;
   sprintf(buffer, "%s.%d", argv[2], myid) ;
   if ( (msgFile = fopen(buffer, "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], argv[2]) ;
      return(-1) ;
   }
   CVfree(buffer) ;
}
type  = atoi(argv[3]) ;
nvtx  = atoi(argv[4]) ;
nitem = atoi(argv[5]) ;
root  = atoi(argv[6]) ;
seed  = atoi(argv[7]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl  -- %d" 
        "\n msgFile -- %s" 
        "\n type    -- %d" 
        "\n nvtx    -- %d" 
        "\n nitem   -- %d" 
        "\n root    -- %d" 
        "\n seed    -- %d" 
        "\n",
        argv[0], msglvl, argv[2], type, nvtx, nitem, root, seed) ;
fflush(msgFile) ;
/*
   -----------------------
   set up the Graph object
   -----------------------
*/
MARKTIME(t1) ;
graph = Graph_new() ;
if ( myid == root ) {
   InpMtx   *inpmtx ;
   int      nedges, totewght, totvwght, v ;
   int      *adj, *vwghts ;
   IVL      *adjIVL, *ewghtIVL ;
/*
   -----------------------
   generate a random graph
   -----------------------
*/
   inpmtx = InpMtx_new() ;
   InpMtx_init(inpmtx, INPMTX_BY_ROWS, INPMTX_INDICES_ONLY, nitem, 0) ;
   Drand_setDefaultFields(&drand) ;
   Drand_setSeed(&drand, seed) ;
   Drand_setUniform(&drand, 0, nvtx) ;
   Drand_fillIvector(&drand, nitem, InpMtx_ivec1(inpmtx)) ;
   Drand_fillIvector(&drand, nitem, InpMtx_ivec2(inpmtx)) ;
   InpMtx_setNent(inpmtx, nitem) ;
   InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n inpmtx mtx filled with raw entries") ;
      InpMtx_writeForHumanEye(inpmtx, msgFile) ;
      fflush(msgFile) ;
   }
   adjIVL = InpMtx_fullAdjacency(inpmtx) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n full adjacency structure") ;
      IVL_writeForHumanEye(adjIVL, msgFile) ;
      fflush(msgFile) ;
   }
   nedges = adjIVL->tsize ;
   if ( type == 1 || type == 3 ) {
      Drand_setUniform(&drand, 1, 10) ;
      vwghts = IVinit(nvtx, 0) ;
      Drand_fillIvector(&drand, nvtx, vwghts) ;
      totvwght = IVsum(nvtx, vwghts) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n vertex weights") ;
         IVfprintf(msgFile, nvtx, vwghts) ;
         fflush(msgFile) ;
      }
   } else {
      vwghts = NULL ;
      totvwght = nvtx ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n totvwght %d", totvwght) ;
      fflush(msgFile) ;
   }
   if ( type == 2 || type == 3 ) {
      ewghtIVL = IVL_new() ;
      IVL_init1(ewghtIVL, IVL_CHUNKED, nvtx) ;
      Drand_setUniform(&drand, 1, 100) ;
      totewght = 0 ;
      for ( v = 0 ; v < nvtx ; v++ ) {
         IVL_listAndSize(adjIVL, v, &size, &adj) ;
         IVL_setList(ewghtIVL, v, size, NULL) ;
         IVL_listAndSize(ewghtIVL, v, &size, &adj) ;
         Drand_fillIvector(&drand, size, adj) ;
         totewght += IVsum(size, adj) ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n ewghtIVL") ;
         IVL_writeForHumanEye(ewghtIVL, msgFile) ;
         fflush(msgFile) ;
      }
   } else {
      ewghtIVL = NULL ;
      totewght = nedges ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n totewght %d", totewght) ;
      fflush(msgFile) ;
   }
   Graph_init2(graph, type, nvtx, 0, nedges, totvwght, totewght,
               adjIVL, vwghts, ewghtIVL) ;
   InpMtx_free(inpmtx) ;
}
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n CPU %8.3f : initialize the Graph object", t2 - t1) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
fflush(msgFile) ;
if ( myid == root ) {
/*
   ----------------------------------------
   compute the checksum of the Graph object
   ----------------------------------------
*/
   chksum = graph->type + graph->nvtx + graph->nvbnd 
          + graph->nedges + graph->totvwght + graph->totewght ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      IVL_listAndSize(graph->adjIVL, v, &size, &list) ;
      chksum += 1 + v + size + IVsum(size, list) ;
   }
   if ( graph->vwghts != NULL ) {
      chksum += IVsum(nvtx, graph->vwghts) ;
   }
   if ( graph->ewghtIVL != NULL ) {
      for ( v = 0 ; v < nvtx ; v++ ) {
         IVL_listAndSize(graph->ewghtIVL, v, &size, &list) ;
         chksum += 1 + v + size + IVsum(size, list) ;
      }
   }
   fprintf(msgFile, "\n\n local chksum = %12.4e", chksum) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   broadcast the Graph object
   --------------------------
*/
MARKTIME(t1) ;
graph = Graph_MPI_Bcast(graph, root, msglvl, msgFile, MPI_COMM_WORLD) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : broadcast the Graph object", t2 - t1) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
/*
   ----------------------------------------
   compute the checksum of the Graph object
   ----------------------------------------
*/
chksum = graph->type + graph->nvtx + graph->nvbnd 
       + graph->nedges + graph->totvwght + graph->totewght ;
for ( v = 0 ; v < nvtx ; v++ ) {
   IVL_listAndSize(graph->adjIVL, v, &size, &list) ;
   chksum += 1 + v + size + IVsum(size, list) ;
}
if ( graph->vwghts != NULL ) {
   chksum += IVsum(nvtx, graph->vwghts) ;
}
if ( graph->ewghtIVL != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      IVL_listAndSize(graph->ewghtIVL, v, &size, &list) ;
      chksum += 1 + v + size + IVsum(size, list) ;
   }
}
fprintf(msgFile, "\n\n local chksum = %12.4e", chksum) ;
fflush(msgFile) ;
/*
   ---------------------------------------
   gather the checksums from the processes
   ---------------------------------------
*/
sums = DVinit(nproc, 0.0) ;
MPI_Gather((void *) &chksum, 1, MPI_DOUBLE, 
           (void *) sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
if ( myid == 0 ) {
   fprintf(msgFile, "\n\n sums") ;
   DVfprintf(msgFile, nproc, sums) ;
   for ( iproc = 0 ; iproc < nproc ; iproc++ ) {
      sums[iproc] -= chksum ;
   }
   fprintf(msgFile, "\n\n errors") ;
   DVfprintf(msgFile, nproc, sums) ;
   fprintf(msgFile, "\n\n maxerror = %12.4e", DVmax(nproc, sums, &loc));
}
/*
   ----------------
   free the objects
   ----------------
*/
DVfree(sums) ;
Graph_free(graph) ;
/*
   ------------------------
   exit the MPI environment
   ------------------------
*/
MPI_Finalize() ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(0) ; }

/*--------------------------------------------------------------------*/
