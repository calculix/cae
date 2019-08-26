/*  Graph_Bcast.c  */

#include "../spoolesMPI.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- to broadcast a Graph IVL object from 
              one processor to all the others

   created -- 98sep10, cca
   -----------------------------------------------
*/
Graph *
Graph_MPI_Bcast (
   Graph      *graph,
   int        root,
   int        msglvl,
   FILE       *msgFile,
   MPI_Comm   comm
) {
int   myid ;
int   itemp[6] ;
/*
   -------------
   find identity
   -------------
*/
MPI_Comm_rank(comm, &myid) ;
if ( myid == root ) {
/*
   -------------------------------
   broadcast the six scalar values
   -------------------------------
*/
   itemp[0] = graph->type     ;
   itemp[1] = graph->nvtx     ;
   itemp[2] = graph->nvbnd    ;
   itemp[3] = graph->nedges   ;
   itemp[4] = graph->totvwght ;
   itemp[5] = graph->totewght ;
   MPI_Bcast((void *) itemp, 6, MPI_INT, root, comm) ;
/*
   ------------------------
   broadcast the IVL object
   ------------------------
*/
   IVL_MPI_Bcast(graph->adjIVL, root, msglvl, msgFile, comm) ;
   if ( graph->type == 1 || graph->type == 3 ) {
/*
      ----------------------------
      broadcast the vertex weights
      ----------------------------
*/
      MPI_Bcast((void *) graph->vwghts, graph->nvtx, MPI_INT, 
                 root, comm) ;
   }
   if ( graph->type == 2 || graph->type == 3 ) {
/*
      -------------------------------------
      broadcast the edge weights IVL object
      -------------------------------------
*/
      IVL_MPI_Bcast(graph->ewghtIVL, root, msglvl, msgFile, comm) ;
   }
} else {
   int   nvtx, type ;
   int   *vwghts ;
   IVL   *adjIVL, *ewghtIVL ;
/*
   --------------
   clear the data
   --------------
*/
   Graph_clearData(graph) ;
/*
   -----------------------------
   receive the six scalar values
   -----------------------------
*/
   MPI_Bcast((void *) itemp, 6, MPI_INT, root, comm) ;
   type = itemp[0] ;
   nvtx = itemp[1] ;
/*
   ----------------------
   receive the IVL object
   ----------------------
*/
   adjIVL = IVL_new() ;
   IVL_MPI_Bcast(adjIVL, root, msglvl, msgFile, comm) ;
   if ( type == 1 || type == 3 ) {
/*
      --------------------------
      receive the vertex weights
      --------------------------
*/
      vwghts = IVinit(nvtx, 0) ;
      MPI_Bcast((void *) vwghts, nvtx, MPI_INT, root, comm) ;
   } else {
      vwghts = NULL ;
   }
   if ( type == 2 || type == 3 ) {
/*
      -----------------------------------
      receive the edge weights IVL object
      -----------------------------------
*/
      ewghtIVL = IVL_new() ;
      IVL_MPI_Bcast(ewghtIVL, root, msglvl, msgFile, comm) ;
   } else {
      ewghtIVL = NULL ;
   }
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
   Graph_init2(graph, type, nvtx, itemp[2], itemp[3], itemp[4],
               itemp[5], adjIVL, vwghts, ewghtIVL) ;
}
return(graph) ; }

/*--------------------------------------------------------------------*/
