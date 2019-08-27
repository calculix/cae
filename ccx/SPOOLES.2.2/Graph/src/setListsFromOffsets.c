/*  setListsFromOffsets.c  */

#include "../Graph.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose -- 

   take an adjacency structure in the (offsets[neqns+1], adjncy[*]) 
   form and load the Graph object. note, pointers to the lists are
   set, no new storage is allocated for the adjacency lists.

   however, during the ordering process each adjacency lists 
   may be shuffled.

   g       -- pointer to Graph object, 
              must be initialized with nvtx = neqns
   neqns   -- # of equations
   offsets -- offsets vector
   adjncy  -- big adjacency vector
      note, the adjacency for list v is found in
            adjncy[offsets[v]:offsets[v+1]-1]
      also note, offsets[] and adjncy[] must be zero based,
      if (offsets,adjncy) come from a harwell-boeing file, they use
      the fortran numbering, so each value must be decremented to
      conform with C's zero based numbering

   created -- 96oct24, cca
   -------------------------------------------------------------------
*/
void
Graph_setListsFromOffsets (
   Graph   *g,
   int     neqns,
   int     offsets[],
   int     adjncy[]
) {
int   v, vsize ;
IVL   *adjIVL ;
/*
   ---------------
   check the input
   ---------------
*/
if (  g == NULL 
   || neqns <= 0
   || offsets == NULL
   || adjncy == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Graph_setListsFromOffsets(%p,%d,%p,%p)"
           "\n bad input\n", g, neqns, offsets, adjncy) ;
   exit(-1) ;
}
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
Graph_init1(g, 0, neqns, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
adjIVL = g->adjIVL ;
/*
   -----------------------------
   set the pointers to the lists
   -----------------------------
*/
for ( v = 0 ; v < neqns ; v++ ) {
   if ( (vsize = offsets[v+1] - offsets[v]) > 0 ) {
      IVL_setPointerToList(adjIVL, v, vsize, &adjncy[offsets[v]]) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
