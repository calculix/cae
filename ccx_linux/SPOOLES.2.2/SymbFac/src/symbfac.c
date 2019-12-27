/*  symbfac.c  */

#include "../SymbFac.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- using one graph, create and return an IVL object 
              that contains a symbolic factorization.

   created -- 97aug29, cca
   -----------------------------------------------------------
*/
IVL *
SymbFac_initFromGraph (
   ETree   *etree,
   Graph   *graph
) {
int    bndweight, count, first, front, ierr, ii, intweight, I, J, 
       last, nfromchildren, nint, nfront, nvtx, size, v, w ;
int    *adj, *bndwghts, *fch, *head, *indices, *keys, *link, *marker, 
       *nodwghts, *sib, *vwghts, *vtxToFront ;
IVL    *frontToVtxIVL ;
Tree   *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || graph == NULL
   || graph->nvtx != nvtx ) {
   fprintf(stderr, "\n fatal error in SymbFac_fromGraph(%p,%p)"
           "\n bad input\n", etree, graph) ;
   if ( etree != NULL ) {
      ETree_writeStats(etree, stderr) ;
   }
   if ( graph != NULL ) {
      Graph_writeStats(graph, stderr) ;
   }
   exit(-1) ;
}
vwghts = graph->vwghts ;
/*
   ----------------------------------------------
   initialize the IVL object to hold the symbolic 
   factorization and set up the work vectors
   ----------------------------------------------
*/
frontToVtxIVL = IVL_new() ;
IVL_init1(frontToVtxIVL, IVL_CHUNKED, nfront) ;
marker     = IVinit(nvtx,   -1) ;
keys       = IVinit(nvtx,   -1) ;
indices    = IVinit(nvtx,   -1) ;
head       = IVinit(nfront, -1) ;
link       = IVinit(nvtx,   -1) ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ;
   head[front] = v ;
}
tree = etree->tree ;
fch  = tree->fch   ;
sib  = tree->sib   ;
/*
   --------------------
   loop over the fronts
   --------------------
*/
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n checking out front %d", J) ;
#endif
   count = 0 ;
/*
   -----------------------------------
   load and mark the internal vertices
   -----------------------------------
*/
   intweight = 0 ;
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      marker[v] = J ;
      indices[count++] = v ;
      intweight += (vwghts != NULL) ? vwghts[v] : 1 ;
   }
   nint = count ;
#if MYDEBUG > 0
   fprintf(stdout, "\n internal : ") ;
   IVfp80(stdout, count, indices, 12, &ierr) ;
#endif
/*
   --------------------------------------------------------
   load vertices from the boundaries of the children fronts
   --------------------------------------------------------
*/
   bndweight = 0 ;
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      IVL_listAndSize(frontToVtxIVL, I, &size, &adj) ;
      for ( ii = size - 1 ; ii >= 0 ; ii-- ) {
         v = adj[ii] ;
         if ( vtxToFront[v] > J ) {
            if ( marker[v] != J ) {
               marker[v] = J ;
               indices[count++] = v ;
               bndweight += (vwghts != NULL) ? vwghts[v] : 1 ;
            }
         } else {
            break ;
         }
      }
   }
   nfromchildren = count - nint ;
#if MYDEBUG > 0
   fprintf(stdout, "\n from children : ") ;
   IVfp80(stdout, count - nint, indices + nint, 17, &ierr) ;
#endif
/*
   ---------------------------------------------------------------
   load unmarked edges that are adjacent to vertices in this front
   ---------------------------------------------------------------
*/
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      Graph_adjAndSize(graph, v, &size, &adj) ; 
      for ( ii = 0 ; ii < size ; ii++ ) {
         w = adj[ii] ;
         if ( w < nvtx && vtxToFront[w] > J && marker[w] != J ) {
            marker[w] = J ;
            indices[count++] = w ;
            bndweight += (vwghts != NULL) ? vwghts[w] : 1 ;
         }
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n from original edges : ") ;
   IVfp80(stdout, count - nint - nfromchildren, 
          indices + nint + nfromchildren, 23, &ierr) ;
#endif
/*
   ---------------------------------
   set the node and boundary weights
   ---------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, 
           "\n front %d, nint %d, count %d, intweight %d, bndweight %d",
           J, nint, count, intweight, bndweight) ;
#endif
   nodwghts[J] = intweight ;
   bndwghts[J] = bndweight ;
/*
   ----------------------------------------------------
   sort the vertices in ascending order of their fronts
   ----------------------------------------------------
*/
   for ( ii = 0 ; ii < count ; ii++ ) {
      v = indices[ii] ;
      keys[ii] = vtxToFront[v] ;
   }
   IV2qsortUp(count, keys, indices) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n indices sorted by fronts") ;
   IVfp80(stdout, count, indices, 25, &ierr) ;
#endif
/*
   -----------------------------------------------------
   sort the indices in ascending order within each front
   -----------------------------------------------------
*/
   front = J ;
   first = 0 ;
   ii    = 1 ;
   while ( ii < count ) {
      v = indices[ii] ;
      if ( vtxToFront[v] != front ) {
         last = ii - 1 ;
         if ( (size = last - first + 1) > 1 ) {
            IVqsortUp(size, indices + first) ;
         }
         first = ii ;
         front = vtxToFront[v] ;
      }
      ii++ ;
   }
   last = count - 1 ;
   if ( (size = last - first + 1) > 1 ) {
      IVqsortUp(size, indices + first) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n indices sorted inside fronts") ;
   IVfp80(stdout, count, indices, 25, &ierr) ;
#endif
/*
   -----------------
   store the indices
   -----------------
*/
      IVL_setList(frontToVtxIVL, J, count, indices) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(indices) ;
IVfree(marker)  ;
IVfree(keys)    ;
IVfree(head)    ;
IVfree(link)    ;

return(frontToVtxIVL) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- create and return an IVL object that 
              contains a symbolic factorization.
   note: we assume that both the ETree and InpMtx objects
         are in their final permutation

   created -- 97mar15, cca
   -------------------------------------------------------
*/
IVL *
SymbFac_initFromInpMtx (
   ETree    *etree,
   InpMtx   *inpmtx
) {
int      count, front, ierr, ii, I, J, nfromchildren, nint, nfront, 
         nvector, nvtx, size, v, w ;
int      *adj, *bndwghts, *fch, *head, *indices, *keys, *link, *marker, 
         *nodwghts, *sib, *vtxToFront ;
IVL      *frontToVtxIVL ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || inpmtx == NULL ) {
   fprintf(stderr, "\n fatal error in Symbfac_initFromInpMtx(%p,%p)"
           "\n bad input\n", etree, inpmtx) ;
   if ( etree != NULL ) {
      ETree_writeStats(etree, stderr) ;
   }
   if ( inpmtx != NULL ) {
      InpMtx_writeStats(inpmtx, stderr) ;
   }
   exit(-1) ;
}
/*
   ------------------------------------------
   check the coordinate type and storage mode
   ------------------------------------------
*/
if ( ! INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   fprintf(stderr, "\n fatal error in Symbfac_initFromInpMtx()"
           "\n bad input, coordType %d, must be INPMTX_BY_CHEVRONS\n",
           InpMtx_coordType(inpmtx)) ;
   exit(-1) ;
}
if ( ! INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(stderr, "\n fatal error in Symbfac_initFromInpMtx()"
           "\n bad input, storageMode %d, must be INPMTX_BY_VECTORS\n",
           InpMtx_storageMode(inpmtx)) ;
   exit(-1) ;
}
nvector = InpMtx_nvector(inpmtx) ;
#if MYDEBUG > 0
fprintf(stdout, "\n nvector = %d", nvector) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------
   initialize the IVL object to hold the symbolic 
   factorization and set up the work vectors
   ----------------------------------------------
*/
frontToVtxIVL = IVL_new() ;
IVL_init1(frontToVtxIVL, IVL_CHUNKED, nfront) ;
marker     = IVinit(nvtx,   -1) ;
keys       = IVinit(nvtx,   -1) ;
indices    = IVinit(nvtx,   -1) ;
head       = IVinit(nfront, -1) ;
link       = IVinit(nvtx,   -1) ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ;
   head[front] = v ;
}
tree = etree->tree ;
fch  = tree->fch   ;
sib  = tree->sib   ;
/*
   --------------------
   loop over the fronts
   --------------------
*/
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n checking out front %d", J) ;
#endif
   count = 0 ;
/*
   -----------------------------------
   load and mark the internal vertices
   -----------------------------------
*/
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      marker[v] = J ;
      indices[count++] = v ;
   }
   nint = count ;
#if MYDEBUG > 0
   fprintf(stdout, "\n internal : ") ;
   IVfp80(stdout, count, indices, 12, &ierr) ;
#endif
/*
   --------------------------------------------------------
   load vertices from the boundaries of the children fronts
   --------------------------------------------------------
*/
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      IVL_listAndSize(frontToVtxIVL, I, &size, &adj) ;
      for ( ii = size - 1 ; ii >= 0 ; ii-- ) {
         v = adj[ii] ;
         if ( vtxToFront[v] > J ) {
            if ( marker[v] != J ) {
               marker[v] = J ;
               indices[count++] = v ;
            }
         } else {
            break ;
         }
      }
   }
   nfromchildren = count - nint ;
#if MYDEBUG > 0
   fprintf(stdout, "\n from children : ") ;
   IVfp80(stdout, count - nint, indices + nint, 17, &ierr) ;
#endif
/*
   ---------------------------------------------------------------
   load unmarked edges that are adjacent to vertices in this front
   ---------------------------------------------------------------
*/
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      if ( v < nvector ) {
         InpMtx_vector(inpmtx, v, &size, &adj) ;
         for ( ii = 0 ; ii < size ; ii++ ) {
            if ( adj[ii] >= 0 ) {
               w = v + adj[ii] ;
            } else {
               w = v - adj[ii] ;
            }
            if ( vtxToFront[w] > J && marker[w] != J ) {
               marker[w] = J ;
               indices[count++] = w ;
            }
         }
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n from original edges : ") ;
   IVfp80(stdout, count - nint - nfromchildren, 
          indices + nint + nfromchildren, 23, &ierr) ;
#endif
/*
   ---------------------------------
   set the node and boundary weights
   ---------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n front %d, nint %d, count %d", J, nint, count) ;
#endif
   nodwghts[J] = nint  ;
   bndwghts[J] = count - nint ;
/*
   ------------------------------------
   sort the vertices in ascending order
   ------------------------------------
*/
   IVqsortUp(count, indices) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n indices sorted ") ;
   IVfp80(stdout, count, indices, 25, &ierr) ;
#endif
/*
   -----------------
   store the indices
   -----------------
*/
   IVL_setList(frontToVtxIVL, J, count, indices) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(indices) ;
IVfree(marker)  ;
IVfree(keys)    ;
IVfree(head)    ;
IVfree(link)    ;

return(frontToVtxIVL) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   compute the symbolic factorization of a matrix pencil
   note: we assume that both the ETree and Pencil 
         objects are in their final permutation

   created -- 98may04, cca
   -----------------------------------------------------
*/
IVL *
SymbFac_initFromPencil (
   ETree    *etree,
   Pencil   *pencil
) {
InpMtx   *inpmtxA, *inpmtxB ;
int      count, front, ierr, ii, I, J, nfromchildren, nint, nfront, 
         nvectorA, nvectorB, nvtx, size, v, w ;
int      *adj, *bndwghts, *fch, *head, *indices, *keys, *link, 
         *marker, *nodwghts, *sib, *vtxToFront ;
IVL      *frontToVtxIVL ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || pencil == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SymbFac_initFromPencil(%p,%p)"
           "\n bad input\n", etree, pencil) ;
   if ( etree != NULL ) {
      ETree_writeStats(etree, stderr) ;
   }
   if ( pencil != NULL ) {
      Pencil_writeStats(pencil, stderr) ;
   }
   exit(-1) ;
}
inpmtxA = pencil->inpmtxA ;
inpmtxB = pencil->inpmtxB ;
/*
   ------------------------------------------
   check the coordinate type and storage mode
   ------------------------------------------
*/
if ( inpmtxA != NULL ) {
   if ( ! INPMTX_IS_BY_CHEVRONS(inpmtxA) ) {
      fprintf(stderr, "\n fatal error in Symbfac_initFromPencil()"
           "\n bad input, coordType %d, must be INPMTX_BY_CHEVRONS\n",
           InpMtx_coordType(inpmtxA)) ;
      exit(-1) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(inpmtxA) ) {
      fprintf(stderr, "\n fatal error in Symbfac_initFromPencil()"
           "\n bad input, storageMode %d, must be INPMTX_BY_VECTORS\n",
           InpMtx_storageMode(inpmtxA)) ;
      exit(-1) ;
   }
   nvectorA = InpMtx_nvector(inpmtxA) ;
} else {
   nvectorA = 0 ;
}
if ( inpmtxB != NULL ) {
   if ( ! INPMTX_IS_BY_CHEVRONS(inpmtxB) ) {
      fprintf(stderr, "\n fatal error in Symbfac_initFromPencil()"
           "\n bad input, coordType %d, must be INPMTX_BY_CHEVRONS\n",
           InpMtx_coordType(inpmtxB)) ;
      exit(-1) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(inpmtxB) ) {
      fprintf(stderr, "\n fatal error in Symbfac_initFromPencil()"
           "\n bad input, storageMode %d, must be INPMTX_BY_VECTORS\n",
           InpMtx_storageMode(inpmtxB)) ;
      exit(-1) ;
   }
   nvectorB = InpMtx_nvector(inpmtxB) ;
} else {
   nvectorB = 0 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n nvectorA = %d", nvectorA) ;
fprintf(stdout, "\n nvectorB = %d", nvectorB) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------
   initialize the IVL object to hold the symbolic 
   factorization and set up the work vectors
   ----------------------------------------------
*/
frontToVtxIVL = IVL_new() ;
IVL_init1(frontToVtxIVL, IVL_CHUNKED, nfront) ;
marker     = IVinit(nvtx,   -1) ;
keys       = IVinit(nvtx,   -1) ;
indices    = IVinit(nvtx,   -1) ;
head       = IVinit(nfront, -1) ;
link       = IVinit(nvtx,   -1) ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   link[v] = head[front] ;
   head[front] = v ;
}
tree = etree->tree ;
fch  = tree->fch   ;
sib  = tree->sib   ;
/*
   --------------------
   loop over the fronts
   --------------------
*/
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n\n checking out front %d", J) ;
#endif
   count = 0 ;
/*
   -----------------------------------
   load and mark the internal vertices
   -----------------------------------
*/
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      marker[v] = J ;
      indices[count++] = v ;
   }
   nint = count ;
#if MYDEBUG > 0
   fprintf(stdout, "\n internal : ") ;
   IVfp80(stdout, count, indices, 12, &ierr) ;
#endif
/*
   --------------------------------------------------------
   load vertices from the boundaries of the children fronts
   --------------------------------------------------------
*/
   for ( I = fch[J] ; I != -1 ; I = sib[I] ) {
      IVL_listAndSize(frontToVtxIVL, I, &size, &adj) ;
      for ( ii = size - 1 ; ii >= 0 ; ii-- ) {
         v = adj[ii] ;
         if ( vtxToFront[v] > J ) {
            if ( marker[v] != J ) {
               marker[v] = J ;
               indices[count++] = v ;
            }
         } else {
            break ;
         }
      }
   }
   nfromchildren = count - nint ;
#if MYDEBUG > 0
   fprintf(stdout, "\n from children : ") ;
   IVfp80(stdout, count - nint, indices + nint, 17, &ierr) ;
#endif
/*
   ---------------------------------------------------------------
   load unmarked edges that are adjacent to vertices in this front
   ---------------------------------------------------------------
*/
   for ( v = head[J] ; v != -1 ; v = link[v] ) {
      if ( inpmtxA != NULL ) {
         InpMtx_vector(inpmtxA, v, &size, &adj) ;
         for ( ii = 0 ; ii < size ; ii++ ) {
            if ( adj[ii] >= 0 ) {
               w = v + adj[ii] ;
            } else {
               w = v - adj[ii] ;
            }
            if ( vtxToFront[w] > J && marker[w] != J ) {
               marker[w] = J ;
               indices[count++] = w ;
            }
         }
      }
      if ( inpmtxB != NULL ) {
         InpMtx_vector(inpmtxB, v, &size, &adj) ;
         for ( ii = 0 ; ii < size ; ii++ ) {
            if ( adj[ii] >= 0 ) {
               w = v + adj[ii] ;
            } else {
               w = v - adj[ii] ;
            }
            if ( vtxToFront[w] > J && marker[w] != J ) {
               marker[w] = J ;
               indices[count++] = w ;
            }
         }
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n from original edges : ") ;
   IVfp80(stdout, count - nint - nfromchildren, 
          indices + nint + nfromchildren, 23, &ierr) ;
#endif
/*
   ---------------------------------
   set the node and boundary weights
   ---------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n front %d, nint %d, count %d", J, nint, count) ;
#endif
   nodwghts[J] = nint  ;
   bndwghts[J] = count - nint ;
/*
   ------------------------------------
   sort the vertices in ascending order
   ------------------------------------
*/
   IVqsortUp(count, indices) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n indices sorted ") ;
   IVfp80(stdout, count, indices, 25, &ierr) ;
#endif
/*
   -----------------
   store the indices
   -----------------
*/
   IVL_setList(frontToVtxIVL, J, count, indices) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(indices) ;
IVfree(marker)  ;
IVfree(keys)    ;
IVfree(head)    ;
IVfree(link)    ;

return(frontToVtxIVL) ; }

/*--------------------------------------------------------------------*/
