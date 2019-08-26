/*  Setup.c  */

#include "../BridgeMPI.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   purpose --

   given an InpMtx object that contains the structure of A, initialize 
     the bridge data structure for the serial factor's and solve's.

   return value --
      1 -- normal return
     -1 -- bridge is NULL
     -2 -- myid is zero and mtxA is NULL

   created -- 98sep25, cca
   -------------------------------------------------------------------
*/
int
BridgeMPI_setup (
   BridgeMPI   *bridge,
   InpMtx      *mtxA
) {
double     t0, t1, t2 ;
FILE       *msgFile ;
int        msglvl, myid, neqns, nproc, rc ;
MPI_Comm   comm ;
/*--------------------------------------------------------------------*/
MARKTIME(t0) ;
/*
   --------------------
   check the input data
   --------------------
*/
if ( bridge == NULL ) {
   fprintf(stderr, "\n fatal error in BridgeMPI_setup()"
           "\n data is NULL\n") ;
   return(-1) ;
}
nproc   = bridge->nproc ;
myid    = bridge->myid  ;
comm    = bridge->comm  ;
msglvl  = bridge->msglvl ;
msgFile = bridge->msgFile ;
neqns   = bridge->neqns ;
if ( myid == 0 ) {
/*
   ------------------------------------------------
   processor 0 does some error checking
   and broadcasts neqns to the other processors
   ------------------------------------------------
*/
   rc = 1 ;
   if ( mtxA == NULL ) {
      fprintf(stderr, "\n fatal error in BridgeMPI_setup()"
              "\n A is NULL\n") ;
      rc = -2 ;
   }
/*
   ------------------------
   broadcast the error flag
   ------------------------
*/
   MPI_Bcast((void *) &rc, 1, MPI_INT, 0, bridge->comm) ;
} else {
/*
   ---------------------------------------
   other processors receive the error flag
   ---------------------------------------
*/
   MPI_Bcast((void *) &rc, 1, MPI_INT, 0, bridge->comm) ;
}
if ( rc != 1 ) {
/*
   ----------------------
   error detected, return
   ----------------------
*/
   return(rc) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   processor 0 does all the non-numeric computation
   (1) orders the graph
   (2) creates the front tree
   (3) creates the symbolic factorization
   ------------------------------------------------
*/
if ( myid == 0 ) {
   ETree    *frontETree ;
   Graph    *graph ;
   int      compressed, nedges, Neqns ;
   IV       *eqmapIV ;
   IVL      *adjIVL, *symbfacIVL ;

   if ( ! (INPMTX_IS_BY_ROWS(mtxA) || INPMTX_IS_BY_COLUMNS(mtxA)) ) {
/*
      ------------------------------
      change coordinate type to rows
      ------------------------------
*/
      InpMtx_changeCoordType(mtxA, INPMTX_BY_ROWS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(mtxA) ) {
/*
      ------------------------------
      change storage mode to vectors
      ------------------------------
*/
      InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
   }
/*
   ---------------------------
   create a Graph object for A
   ---------------------------
*/
   MARKTIME(t1) ;
   graph  = Graph_new() ;
   adjIVL = InpMtx_fullAdjacency(mtxA);
   nedges = bridge->nedges = IVL_tsize(adjIVL),
   Graph_init2(graph, 0, neqns, 0, nedges,
               neqns, nedges, adjIVL, NULL, NULL) ;
   MARKTIME(t2) ;
   bridge->cpus[0] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : time to create Graph", t2 - t1) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n graph of the input matrix") ;
      Graph_writeForHumanEye(graph, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------------
   get the equivalence map
   -----------------------
*/
   MARKTIME(t1) ;
   eqmapIV = Graph_equivMap(graph) ;
   Neqns = bridge->Neqns = 1 + IV_max(eqmapIV) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n graph's equivalence map") ;
      IV_writeForHumanEye(eqmapIV, msgFile) ;
      fflush(msgFile) ;
   }
   if ( Neqns < bridge->compressCutoff * neqns ) {
      Graph   *cgraph ;
/*
      ------------------
      compress the graph
      ------------------
*/
      cgraph = Graph_compress2(graph, eqmapIV, 1) ;
      Graph_free(graph) ;
      graph = cgraph ;
      compressed = 1 ;
      bridge->Nedges = graph->nedges ;
   } else {
      compressed = 0 ;
   }
   MARKTIME(t2) ;
   bridge->cpus[1] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : time to create compressed graph",
              t2 - t1) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n graph to order") ;
      Graph_writeForHumanEye(graph, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ---------------
   order the graph
   ---------------
*/
   MARKTIME(t1) ;
   if ( bridge->maxdomainsize <= 0 ) {
      bridge->maxdomainsize = neqns/32 ;
   }
   if ( bridge->maxdomainsize <= 0 ) {
      bridge->maxdomainsize = 1 ;
   }
   if ( bridge->maxnzeros < 0 ) {
      bridge->maxnzeros = 0.01*neqns ;
   }
   if ( bridge->maxsize < 0 ) {
      bridge->maxsize = 64 ;
   }
   frontETree = orderViaBestOfNDandMS(graph, bridge->maxdomainsize,
                                     bridge->maxnzeros, bridge->maxsize,
                                     bridge->seed, msglvl, msgFile) ;
   MARKTIME(t2) ;
   bridge->cpus[2] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : time to order graph", t2 - t1) ;
      fflush(msgFile) ;
   }
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n front tree from ordering") ;
      ETree_writeForHumanEye(frontETree, msgFile) ;
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   if ( compressed == 1 ) {
      ETree   *etree ;
      IVL     *tempIVL ;
/*
      ----------------------------------------------------------
      compute the symbolic factorization of the compressed graph
      ----------------------------------------------------------
*/
      tempIVL = SymbFac_initFromGraph(frontETree, graph) ;
/*
      -------------------------------------------------------
      expand the symbolic factorization to the original graph
      -------------------------------------------------------
*/
      symbfacIVL = IVL_expand(tempIVL, eqmapIV) ;
      IVL_free(tempIVL) ;
/*
      ---------------------
      expand the front tree
      ---------------------
*/
      etree = ETree_expand(frontETree, eqmapIV) ;
      ETree_free(frontETree) ;
      frontETree = etree ;
   } else {
/*
      --------------------------------------------------------
      compute the symbolic factorization of the original graph
      --------------------------------------------------------
*/
      symbfacIVL = SymbFac_initFromGraph(frontETree, graph) ;
   }
   MARKTIME(t2) ;
   bridge->frontETree = frontETree ;
   bridge->symbfacIVL = symbfacIVL ;
/*
   ----------------------------------------------
   get the old-to-new and new-to-old permutations
   ----------------------------------------------
*/
   bridge->oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
   bridge->newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n old-to-new permutation") ;
      IV_writeForHumanEye(bridge->oldToNewIV, msgFile) ;
      fprintf(msgFile, "\n\n new-to-old permutation") ;
      IV_writeForHumanEye(bridge->newToOldIV, msgFile) ;
      fflush(msgFile) ;
   }
/*
   ------------------------------------------------------
   overwrite the symbolic factorization with the permuted 
   indices and sort the lists into ascending order
   ------------------------------------------------------
*/
   IVL_overwrite(symbfacIVL, bridge->oldToNewIV) ;
   IVL_sortUp(symbfacIVL) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n symbolic factorization") ;
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
      fflush(msgFile) ;
   }
/*
   --------------------------------------
   permute the vertices in the front tree
   --------------------------------------
*/
   ETree_permuteVertices(frontETree, bridge->oldToNewIV) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n permuted front etree") ;
      ETree_writeForHumanEye(frontETree, msgFile) ;
      fflush(msgFile) ;
   }
   MARKTIME(t2) ;
   bridge->cpus[3] += t2 - t1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n CPU %8.3f : time for symbolic factorization",
              t2 - t1) ;
      fflush(msgFile) ;
   }
/*
   ------------------------
   free the working storage
   ------------------------
*/
   Graph_free(graph) ;
   IV_free(eqmapIV) ;
}

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   processor 0 broadcasts the ETree and the symbolic factorization IVL 
   -------------------------------------------------------------------
*/
MPI_Barrier(comm) ;
if ( myid == 0 ) {
   MARKTIME(t1) ;
   ETree_MPI_Bcast(bridge->frontETree, 0, msglvl, msgFile, comm) ;
   MARKTIME(t2) ;
   bridge->cpus[4] += t2 - t1 ;
   t1 = t2 ;
   IVL_MPI_Bcast(bridge->symbfacIVL, 0, msglvl, msgFile, comm) ;
   MARKTIME(t2) ;
   bridge->cpus[5] += t2 - t1 ;
} else {
   MARKTIME(t1) ;
   bridge->frontETree = ETree_MPI_Bcast(NULL, 0, msglvl, msgFile, comm);
   MARKTIME(t2) ;
   bridge->cpus[4] += t2 - t1 ;
   t1 = t2 ;
   bridge->symbfacIVL = IVL_MPI_Bcast(NULL, 0, msglvl, msgFile, comm) ;
   MARKTIME(t2) ;
   bridge->cpus[5] += t2 - t1 ;
}
MARKTIME(t2) ;
bridge->cpus[6] += t2 - t0 ;
/*--------------------------------------------------------------------*/

return(1) ; }

/*--------------------------------------------------------------------*/
