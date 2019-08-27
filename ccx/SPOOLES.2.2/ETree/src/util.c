/*  util.c  */

#include "../ETree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95nov15, cca
   ----------------------------------------------
*/
int
ETree_sizeOf (
   ETree   *etree
) {
int   nbytes ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_sizeOf(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
nbytes = sizeof(struct _ETree) ;
if ( etree->tree != NULL ) {
   nbytes += Tree_sizeOf(etree->tree) ;
}
if ( etree->nodwghtsIV != NULL ) {
   nbytes += IV_sizeOf(etree->nodwghtsIV) ;
}
if ( etree->nodwghtsIV != NULL ) {
   nbytes += IV_sizeOf(etree->bndwghtsIV) ;
}
if ( etree->vtxToFrontIV != NULL ) {
   nbytes += IV_sizeOf(etree->vtxToFrontIV) ;
}
return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   return the number of factor indices

   created  -- 95nov15, cca
   modified -- 96jan11, cca
   ----------------------------------------
*/
int
ETree_nFactorIndices (
   ETree   *etree
) {
int   nb, nfront, nv, nind, nvtx, v ;
int   *bndwghts, *nodwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nFactorIndices(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
nind = 0 ;
for ( v = 0 ; v < nfront ; v++ ) {
   nv = nodwghts[v] ;
   nb = bndwghts[v] ;
   nind += nv + nb ;
}
return(nind) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return the number of factor entries

   symflag -- symmetry flag
     0 (SPOOLES_SYMMETRIC)    -- symmetric
     1 (SPOOLES_HERMITIAN)    -- hermitian
     2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric

   created -- 98jun05, cca
   ------------------------------------------
*/
int
ETree_nFactorEntries (
   ETree   *etree,
   int     symflag
) {
int   J, nfront, nvtx, nzf ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nFactorEntries(%p,%d)"
           "\n bad input\n", etree, symflag) ;
   exit(-1) ;
}
nzf = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   nzf += ETree_nFactorEntriesInFront(etree, symflag, J) ;
}
return(nzf) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return the number of factor operations

   type -- type of matrix entries
     1 (SPOOLES_REAL)    -- real entries
     2 (SPOOLES_COMPLEX) -- complex entries
   symflag -- symmetry flag
     0 (SPOOLES_SYMMETRIC)    -- symmetric
     1 (SPOOLES_HERMITIAN)    -- hermitian
     2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric

   created -- 98jun05, cca
   ------------------------------------------
*/
double
ETree_nFactorOps (
   ETree   *etree,
   int     type,
   int     symflag
) {
double   ops ;
int      J, nfront, nvtx ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nFactorOps(%p,%d,%d)"
           "\n bad input\n", etree, type, symflag) ;
   exit(-1) ;
}
ops = 0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   ops += ETree_nInternalOpsInFront(etree, type, symflag, J)
        + ETree_nExternalOpsInFront(etree, type, symflag, J) ;
}
return(ops) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   return the number of entries an LU front

   created -- 96dec04, cca
   ----------------------------------------
*/
double
ETree_nFactorEntriesInFront (
   ETree   *etree,
   int     symflag, 
   int     J 
) {
int   b, m, nent ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || etree->nfront <= 0
   || J < 0 || J >= etree->nfront ) {
   fprintf(stderr, 
           "\n fatal error in ETree_nFactorEntriesInFront(%p,%d,%d)"
           "\n bad input\n", etree, symflag, J) ;
   exit(-1) ;
}
b = IV_entry(etree->nodwghtsIV, J) ;
m = IV_entry(etree->bndwghtsIV, J) ;
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
   nent = (b*(b+1))/2 + b*m ;
   break ;
case SPOOLES_NONSYMMETRIC :
   nent = b*b + 2*b*m ;
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in ETree_nFactorEntriesInFront(%p,%d,%d)"
           "\n bad symflag\n", etree, symflag, J) ;
   break ;
}

return(nent) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   return the number of internal LU operations for a front

   created -- 96dec04, cca
   -------------------------------------------------------
*/
double
ETree_nInternalOpsInFront (
   ETree   *etree,
   int     type, 
   int     symflag, 
   int     J 
) {
double   b, m, ops ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || etree->nfront <= 0
   || J < 0 || J >= etree->nfront ) {
   fprintf(stderr, 
           "\n fatal error in ETree_nInternalOpsInFront(%p,%d,%d,%d)"
           "\n bad input\n", etree, type, symflag, J) ;
   exit(-1) ;
}
b = ETree_frontSize(etree, J) ;
m = ETree_frontBoundarySize(etree, J) ;
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
   ops = (b*(b+1)*(2*b+1))/6. + m*b*b ;
   break ;
case SPOOLES_NONSYMMETRIC :
   ops = b*(2*b*b+1)/3. + 2*m*b*b ;
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in ETree_nInternalOpsInFront(%p,%d,%d,%d)"
           "\n bad symflag\n", etree, type, symflag, J) ;
   break ;
}
switch ( type ) {
case SPOOLES_REAL :
   break ;
case SPOOLES_COMPLEX :
   ops = 4*ops ;
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in ETree_nInternalOpsInFront(%p,%d,%d,%d)"
           "\n bad type\n", etree, type, symflag, J) ;
   break ;
}
return(ops) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   return the number of external LU operations for a front

   created -- 96dec04, cca
   -------------------------------------------------------
*/
double
ETree_nExternalOpsInFront (
   ETree   *etree,
   int     type,
   int     symflag,
   int     J 
) {
double   b, m, ops ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL 
   || etree->nfront <= 0
   || J < 0 || J >= etree->nfront ) {
   fprintf(stderr, 
           "\n fatal error in ETree_nExternalOpsInFront(%p,%d,%d,%d)"
           "\n bad input\n", etree, J, type, symflag) ;
   exit(-1) ;
}
b = IV_entry(etree->nodwghtsIV, J) ;
m = IV_entry(etree->bndwghtsIV, J) ;
switch ( symflag ) {
case SPOOLES_SYMMETRIC :
case SPOOLES_HERMITIAN :
   ops = m*(m+1)*b ;
   break ;
case SPOOLES_NONSYMMETRIC :
   ops = 2*b*m*m ;
   break ;
default :
   break ;
}
switch ( type ) {
case SPOOLES_REAL :
   break ;
case SPOOLES_COMPLEX :
   ops = 4*ops ;
   break ;
default :
   fprintf(stderr, 
           "\n fatal error in ETree_nExternalOpsInFront(%p,%d,%d,%d)"
           "\n bad input\n", etree, J, type, symflag) ;
   break ;
}
return(ops) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   return a DV object that contains the number of operations 
   for each front using a backward looking algorithm

   created -- 96dec04, cca
   ---------------------------------------------------------
*/
DV *
ETree_backwardOps (
   ETree   *etree,
   int     type,
   int     symflag,
   int     *vwghts,
   IVL     *symbfacIVL
) {
double   extops, opsKbndK, opsKK ;
double   *ops ;
DV       *opsDV ;
int      bndwghtJ, ii, J, k, K, kwght, 
         nadj, nfront, size, wghtJ, wghtK ;
int      *counts, *indices, *list, *mark, *vtxToFront ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || symbfacIVL == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_backwardOps(%p,%p,%p)"
           "\n bad input\n", etree, vwghts, symbfacIVL) ;
   exit(-1) ;
}
nfront     = etree->nfront ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
list   = IVinit(nfront, -1) ;
mark   = IVinit(nfront, -1) ;
counts = IVinit(nfront, 0) ;
/*
   ----------------------------
   initialize the ops DV object
   ----------------------------
*/
opsDV = DV_new() ;
DV_init(opsDV, nfront, NULL) ;
ops = DV_entries(opsDV) ;
DV_fill(opsDV, 0.0) ;
/*
   -------------------
   fill the ops vector
   -------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   ops[J]  += ETree_nInternalOpsInFront(etree, type, symflag, J) ;
   wghtJ    = ETree_frontSize(etree, J) ;
   bndwghtJ = ETree_frontBoundarySize(etree, J) ;
   IVL_listAndSize(symbfacIVL, J, &size, &indices) ;
   for ( ii = nadj = 0 ; ii < size ; ii++ ) {
      k = indices[ii] ;
      if ( (K = vtxToFront[k]) != J ) {
         kwght = (vwghts == NULL) ? 1 : vwghts[k] ;
         if ( mark[K] != J ) {
            counts[K] = 0 ;
            mark[K] = J ;
            list[nadj++] = K ;
         }
         counts[K] += kwght ;
      }
   }
   IVqsortUp(nadj, list) ;
   extops = 0.0 ;
   for ( ii = 0 ; ii < nadj ; ii++ ) {
      K = list[ii] ;
      wghtK = counts[K] ;
      bndwghtJ -= wghtK ;
      if ( type == SPOOLES_REAL ) {
         opsKbndK = 2*wghtJ*wghtK*bndwghtJ ;
         if ( symflag == SPOOLES_SYMMETRIC ) {
            opsKK = wghtJ*wghtK*(wghtK+1) ;
         } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
            opsKK = 2*wghtJ*wghtK*wghtK ;
         }
      } else if ( type == SPOOLES_COMPLEX ) {
         opsKbndK = 8*wghtJ*wghtK*bndwghtJ ;
         if (  symflag == SPOOLES_SYMMETRIC 
            || symflag == SPOOLES_HERMITIAN ) {
            opsKK = 4*wghtJ*wghtK*(wghtK+1) ;
         } else if ( symflag == SPOOLES_NONSYMMETRIC ) {
            opsKK = 8*wghtJ*wghtK*wghtK ;
         }
      }
      extops += opsKK + opsKbndK ;
      ops[K] += opsKK + opsKbndK ;
      if ( symflag == SPOOLES_NONSYMMETRIC ) {
         extops += opsKbndK ;
         ops[K] += opsKbndK ;
      }
   }
/*
   fprintf(stdout, 
       "\n front %d, %.0f internal ops, %.0f external ops, %.0f extops",
           J, ETree_nInternalOpsInFront(etree, type, symflag, J),
           ETree_nExternalOpsInFront(etree, type, symflag, J), extops) ;
*/
}
IVfree(list) ;
IVfree(mark) ;
IVfree(counts) ;

return(opsDV) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   return an IV object that contains 
   the number of entries for each front 

   created -- 98jan30, cca
   ------------------------------------
*/
IV *
ETree_factorEntriesIV (
   ETree   *etree,
   int     symflag
) {
int   J, nfront ;
int   *nzf ;
IV    *nzfIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_factorEntriesIV(%p,%d)"
           "\n bad input\n", etree, symflag) ;
   exit(-1) ;
}
nfront = ETree_nfront(etree) ;
/*
   ------------------------
   initialize the IV object
   ------------------------
*/
nzfIV = IV_new() ;
IV_init(nzfIV, nfront, NULL) ;
nzf = IV_entries(nzfIV) ;
IV_fill(nzfIV, 0) ;
/*
   -------------------
   fill the nzf vector
   -------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   nzf[J] = ETree_nFactorEntriesInFront(etree, symflag, J) ;
}
return(nzfIV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   return a DV object that contains the number of operations 
   for each front using a forward-looking algorithm

   created -- 96dec04, cca
   ---------------------------------------------------------
*/
DV *
ETree_forwardOps (
   ETree   *etree,
   int     type,
   int     symflag
) {
double   *ops ;
DV       *opsDV ;
int      J, nfront ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_forwardOps(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
nfront     = etree->nfront ;
opsDV = DV_new() ;
DV_init(opsDV, nfront, NULL) ;
ops = DV_entries(opsDV) ;
DV_fill(opsDV, 0.0) ;
for ( J = 0 ; J < nfront ; J++ ) {
   ops[J] += ETree_nInternalOpsInFront(etree, type, symflag, J)
          +  ETree_nExternalOpsInFront(etree, type, symflag, J) ;
}
return(opsDV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   given an IV object that maps uncompressed vertices to vertices,
   create and return an ETree object that is relative to the 
   uncompressed graph.

   created -- 97feb13, cca
   ---------------------------------------------------------------
*/
ETree *
ETree_expand (
   ETree   *etree,
   IV      *eqmapIV
) {
ETree   *etree2 ;
int     ii, ndof, nfront ;
int     *map, *vtxToFront, *vtxToFront2 ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || eqmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_expand(%p,%p)"
           "\n bad input\n", etree, eqmapIV) ;
   exit(-1) ;
}
nfront = etree->nfront ;
IV_sizeAndEntries(eqmapIV, &ndof, &map) ;
/*
   ---------------------------
   create the new ETree object
   ---------------------------
*/
etree2 = ETree_new() ;
ETree_init1(etree2, nfront, ndof) ;
IV_copy(etree2->nodwghtsIV, etree->nodwghtsIV) ;
IV_copy(etree2->bndwghtsIV, etree->bndwghtsIV) ;
etree2->tree->root = etree->tree->root ;
IVcopy(nfront, etree2->tree->par, etree->tree->par) ;
IVcopy(nfront, etree2->tree->fch, etree->tree->fch) ;
IVcopy(nfront, etree2->tree->sib, etree->tree->sib) ;
vtxToFront  = IV_entries(etree->vtxToFrontIV) ;
vtxToFront2 = IV_entries(etree2->vtxToFrontIV) ;
for ( ii = 0 ; ii < ndof ; ii++ ) {
   vtxToFront2[ii] = vtxToFront[map[ii]] ;
}

return(etree2) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   this method is used to splice together two front trees
   when the domain vertices and schur complement vertices
   have been ordered separately.

   etree0 -- the lower front tree is for vertices in the domains.
   graph0 -- graph for all the vertices
   mapIV  -- IV object that maps vertices to schur complement
             vertices, if IV_entry(mapIV, v) < 0 then v is 
             a domain vertex.
   etree1 -- the upper front tree is for vertices in the schur 
             complement.

   created -- 97feb01, cca
   --------------------------------------------------------------
*/
ETree *
ETree_spliceTwoETrees (
   ETree   *etree0,
   Graph   *graph0,
   IV      *mapIV,
   ETree   *etree1
) {
ETree   *etree2 ;
int     *bndwghts0, *bndwghts1, *bndwghts2, *fch0, *head0, *link0, 
        *map, *mark, *nodwghts0, *nodwghts1, *nodwghts2, *par0, 
        *par1, *par2, *sib0, *vadj, *vtxToFront0, *vtxToFront1, 
        *vtxToFront2 ;
int     ii, J, K, nfront0, nfront1, nfront2, nvtx, phi, v, vsize, w ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree0 == NULL || graph0 == NULL 
   || mapIV == NULL || etree1 == NULL ) {
   fprintf(stderr, 
           "\n fatal error in ETree_spliceTwoETrees(%p,%p,%p,%p)"
           "\n bad input\n",
           etree0, graph0, mapIV, etree1) ;
   exit(-1) ;
}
nfront0     = etree0->nfront ;
nvtx        = etree0->nvtx    ;
par0        = etree0->tree->par ;
fch0        = etree0->tree->fch ;
sib0        = etree0->tree->sib ;
nodwghts0   = IV_entries(etree0->nodwghtsIV) ;
bndwghts0   = IV_entries(etree0->bndwghtsIV) ;
vtxToFront0 = IV_entries(etree0->vtxToFrontIV) ;
nfront1     = etree1->nfront ;
par1        = etree1->tree->par ;
bndwghts1   = IV_entries(etree1->bndwghtsIV) ;
nodwghts1   = IV_entries(etree1->nodwghtsIV) ;
vtxToFront1 = IV_entries(etree1->vtxToFrontIV) ;
map         = IV_entries(mapIV) ;
/*
   -------------------------
   create the new front tree
   -------------------------
*/
nfront2 = nfront0 + nfront1 ;
etree2 = ETree_new() ;
ETree_init1(etree2, nfront2, etree0->nvtx) ;
par2        = etree2->tree->par ;
nodwghts2   = IV_entries(etree2->nodwghtsIV) ;
bndwghts2   = IV_entries(etree2->bndwghtsIV) ;
vtxToFront2 = IV_entries(etree2->vtxToFrontIV) ;
/*
   --------------------------------------------------
   fill the parent fields for fronts in the same tree
   --------------------------------------------------
*/
for ( J = 0 ; J < nfront0 ; J++ ) {
   par2[J]      = par0[J] ;
   nodwghts2[J] = nodwghts0[J] ;
   bndwghts2[J] = bndwghts0[J] ;
}
for ( J = 0 ; J < nfront1 ; J++ ) {
   par2[J+nfront0]      = nfront0 + par1[J] ;
   nodwghts2[J+nfront0] = nodwghts1[J] ;
   bndwghts2[J+nfront0] = bndwghts1[J] ;
}
/*
   ---------------------------
   set the vertex to front map
   ---------------------------
*/
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( (J = vtxToFront0[v]) >= 0 ) {
      vtxToFront2[v] = J ;
   } else {
      vtxToFront2[v] = vtxToFront1[map[v]] + nfront0 ;
   }
}
/*
   ---------------------------------------------
   link the vertices to fronts in the lower tree
   ---------------------------------------------
*/
head0 = IVinit(nfront0, -1) ;
link0 = IVinit(nvtx,    -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( (J = vtxToFront0[v]) >= 0 ) {
      link0[v] = head0[J] ;
      head0[J] = v ;
   }
}
/*
   -------------------------------------------------------
   link roots of the lower tree to nodes in the upper tree
   -------------------------------------------------------
*/
mark = IVinit(nvtx, -1) ;
for ( J = etree0->tree->root ; J != -1 ; J = sib0[J] ) {
/*
   ---------------------------------------
   K is the parent front in the upper tree
   ---------------------------------------
*/
   K = nfront1 ;
/*
   ---------------------------------------
   loop over vertices in the lower front J
   ---------------------------------------
*/
   for ( v = head0[J] ; v != -1 ; v = link0[v] ) {
      Graph_adjAndSize(graph0, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( vtxToFront0[w] < 0 ) {
            phi = map[w] ;
/*
            ---------------------------------------------------
            w is a vertex that belongs to phi in the upper tree
            ---------------------------------------------------
*/
            if ( mark[phi] != J ) {
               mark[phi] = J ;
               if ( K > vtxToFront1[phi] ) {
/*
                  ------------------------------------
                  least numbered adjacent front so far
                  ------------------------------------
*/
                  K = vtxToFront1[phi] ;
               }
            }
         }
      }
   }
   if ( K < nfront1 ) {
   /*
      --------------------
      set the parent field
      --------------------
   */
      par2[J] = nfront0 + K ;
   }
}
/*
   -----------------------------
   set the remaining tree fields
   -----------------------------
*/
Tree_setFchSibRoot(etree2->tree) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head0) ;
IVfree(link0) ;
IVfree(mark)  ;

return(etree2) ; }

/*--------------------------------------------------------------------*/
