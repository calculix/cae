/*  ms.c  */

#include "../ETree.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   returns a compidsIV IV object that maps the
   vertices to a domain (compids[v] > 1)
   or to the multisector (compids[v] = 0).
   the vertices in the multisector is specified 
   by their depth of their front in the tree.

   created -- 96jan04, cca
   ------------------------------------------------
*/
IV *
ETree_msByDepth (
   ETree   *etree,
   int     depth
) {
int    front, nfront, nvtx, v ;
int    *compids, *dmetric, *vtxToFront ;
IV     *compidsIV, *vmetricIV, *dmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (   etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0
   || depth <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_msByDepth(%p,%d)"
           "\n bad input\n", etree, depth) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   --------------------
   get the depth metric
   --------------------
*/
vmetricIV = IV_new() ;
IV_init(vmetricIV, nfront, NULL) ;
IV_fill(vmetricIV, 1) ;
dmetricIV = Tree_setDepthImetric(etree->tree, vmetricIV) ;
dmetric   = IV_entries(dmetricIV) ;
#if MYDEBUG > 0
{ int   ierr ;
fprintf(stdout, "\n ETree_msByDepth") ;
fprintf(stdout, "\n vmetric") ;
IV_writeForHumanEye(vmetricIV, stdout) ;
fprintf(stdout, "\n dmetric") ;
IV_writeForHumanEye(dmetricIV, stdout) ;
}
#endif
IV_free(vmetricIV) ;
/*
   --------------------------------------------------
   fill the IV object with the ids in the multisector
   --------------------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, nvtx, NULL) ;
compids = IV_entries(compidsIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   if ( dmetric[front] <= depth ) {
      compids[v] = 0 ;
   } else {
      compids[v] = 1 ;
   }
}
IV_free(dmetricIV) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   construct a multisector based on vertices found in a subtree.

   created -- 96jan04, cca
   ----------------------------------------------------------------
*/
IV *
ETree_msByNvtxCutoff (
   ETree    *etree,
   double   cutoff
) {
int      front, nfront, nvtx, v ;
int      *compids, *tmetric, *vtxToFront ;
IV       *compidsIV, *tmetricIV, *vmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_msByCutoff(%p,%f)"
           "\n bad input\n", etree, cutoff) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------
   get the subtree metric
   ----------------------
*/
vmetricIV = ETree_nvtxMetric(etree) ;
tmetricIV = Tree_setSubtreeImetric(etree->tree, vmetricIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ETree_msByNvtxCutoff") ;
fprintf(stdout, "\n vmetric") ;
IV_writeForHumanEye(vmetricIV, stdout) ;
fprintf(stdout, "\n tmetric") ;
IV_writeForHumanEye(tmetricIV, stdout) ;
fflush(stdout) ;
#endif
IV_free(vmetricIV) ;
cutoff = cutoff * IV_max(tmetricIV) ;
/*
   --------------------------------------------------
   fill the IV object with the ids in the multisector
   --------------------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, nvtx, NULL) ;
compids = IV_entries(compidsIV) ;
tmetric = IV_entries(tmetricIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   if ( tmetric[front] >= cutoff ) {
      compids[v] = 0 ;
   } else {
      compids[v] = 1 ;
   }
}
IV_free(tmetricIV) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   construct a multisector based on the number 
   of factor entries found in a subtree.

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 96jan04, cca
   --------------------------------------------------
*/
IV *
ETree_msByNentCutoff (
   ETree    *etree,
   double   cutoff,
   int      symflag
) {
int      front, nfront, nvtx, v ;
int      *compids, *tmetric, *vtxToFront ;
IV       *compidsIV, *tmetricIV, *vmetricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_msByCutoff(%p,%f,%d)"
           "\n bad input\n", etree, cutoff, symflag) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------
   get the subtree metric
   ----------------------
*/
vmetricIV = ETree_nentMetric(etree, symflag) ;
tmetricIV = Tree_setSubtreeImetric(etree->tree, vmetricIV) ;
#if MYDEBUG > 0
fprintf(stdout, "\n ETree_msByNentCutoff") ;
fprintf(stdout, "\n vmetric") ;
IV_writeForHumanEye(vmetricIV, stdout) ;
fprintf(stdout, "\n tmetric") ;
IV_writeForHumanEye(tmetricIV, stdout) ;
fflush(stdout) ;
#endif
IV_free(vmetricIV) ;
cutoff = cutoff * IV_max(tmetricIV) ;
/*
   --------------------------------------------------
   fill the IV object with the ids in the multisector
   --------------------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, nvtx, NULL) ;
compids = IV_entries(compidsIV) ;
tmetric = IV_entries(tmetricIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   if ( tmetric[front] >= cutoff ) {
      compids[v] = 0 ;
   } else {
      compids[v] = 1 ;
   }
}
IV_free(tmetricIV) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   construct a multisector based on the number 
   of factor operations found in a subtree.

   type -- type of entries,
     SPOOLES_REAL or SPOOLES_COMPLEX

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 96jan04, cca
   --------------------------------------------------
*/
IV *
ETree_msByNopsCutoff (
   ETree    *etree,
   double   cutoff,
   int      type,
   int      symflag
) {
double   *tmetric ;
DV       *tmetricDV, *vmetricDV ;
int      front, nfront, nvtx, v ;
int      *compids, *vtxToFront ;
IV       *compidsIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0
   || (nvtx = etree->nvtx) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_msByCutoff(%p,%f,%d)"
           "\n bad input\n", etree, cutoff, symflag) ;
   exit(-1) ;
}
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------
   get the subtree metric
   ----------------------
*/
vmetricDV = ETree_nopsMetric(etree, type, symflag) ;
tmetricDV = Tree_setSubtreeDmetric(etree->tree, vmetricDV) ;
#if MYDEBUG >= 0
fprintf(stdout, "\n ETree_msByNopsCutoff") ;
fprintf(stdout, "\n vmetric") ;
DV_writeForHumanEye(vmetricDV, stdout) ;
fprintf(stdout, "\n tmetric") ;
DV_writeForHumanEye(tmetricDV, stdout) ;
fflush(stdout) ;
fprintf(stdout, "\n max(tmetricDV) = %.0f, sum(vmetricDV) = %.0f",
        DV_max(tmetricDV), DV_sum(vmetricDV)) ;
fprintf(stdout, "\n cutoff = %.0f", cutoff * DV_max(tmetricDV)) ;
#endif
cutoff = cutoff * DV_max(tmetricDV) ;
/*
   --------------------------------------------------
   fill the IV object with the ids in the multisector
   --------------------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, nvtx, NULL) ;
compids = IV_entries(compidsIV) ;
tmetric = DV_entries(tmetricDV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   front = vtxToFront[v] ;
   if ( tmetric[front] >= cutoff ) {
      compids[v] = 0 ;
   } else {
      compids[v] = 1 ;
   }
}
{
double   domops, schurops ;
double   *vmetric ;
int      J ;

vmetric = DV_entries(vmetricDV) ;
domops = schurops = 0.0 ;
for ( J = 0 ; J < nfront ; J++ ) {
   if ( tmetric[J] >= cutoff ) {
      schurops += vmetric[J] ;
   } else {
      domops += vmetric[J] ;
   }
}
fprintf(stdout, "\n domops = %.0f, schurops = %.0f, total = %.0f",
        domops, schurops, domops + schurops) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(vmetricDV) ;
DV_free(tmetricDV) ;

return(compidsIV) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   purpose -- given a front tree and a multisector map vector,
     fill the map vector with domain ids and the three statistics
     arrays with domain and schur complement statistics.

   frontETree -- front tree object, unchanged on output
   msIV -- map from fronts to domains or schur complement
     on input, ms[v] = 0 --> v is in the schur complement
               ms[v] = 1 --> v is not in the schur complement
     on output, ms[v] =  0 --> v is in the schur complement
                ms[v] != 0 --> v is in domain ms[v]
   on output
      nvtxIV -- nvtx[ireg] = # of dof in region ireg
      nzfIV  -- nzf[ireg] = # of factor entries in region ireg
      opsIV  -- ops[ireg] = # of factor ops in region ireg

   type -- type of entries, SPOOLES_REAL or SPOOLES_COMPLEX

   symflag -- symmetry flag, one of SPOOLES_SYMMETRIC
     SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 98jan30, cca
   --------------------------------------------------------------
*/
void
ETree_msStats (
   ETree   *frontETree,
   IV      *msIV,
   IV      *nvtxIV,
   IV      *nzfIV,
   DV      *opsDV,
   int     type,
   int     symflag
) {
double   *opsreg, *opsvec ;
DV       *opsvecDV ;
int      J, K, ndom, nfront, nvtx, reg, v ;
int      *map, *ms, *nodwghts, *nvtxreg, 
         *nzfreg, *nzfvec, *par, *vtxToFront ;
IV       *nzfvecIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontETree == NULL || msIV == NULL || nvtxIV == NULL
   || nzfIV == NULL || opsDV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_msStats()"
           "\n frontETree = %p, msIV = %p, nvtxIV = %p"
           "\n nzfIV = %p, opsDV = %p, symflag = %d\n",
           frontETree, msIV, nvtxIV, nzfIV, opsDV, symflag) ;
   exit(-1) ;
}
nfront     = ETree_nfront(frontETree) ;
nvtx       = ETree_nvtx(frontETree) ;
tree       = ETree_tree(frontETree) ;
par        = ETree_par(frontETree) ;
vtxToFront = ETree_vtxToFront(frontETree) ;
ms         = IV_entries(msIV) ;
/*
   ----------------------------------------
   if J is not in the schur complement then
      fill ms[J] with its domain id
   endif
   ----------------------------------------
*/
map = IVinit(nfront, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   map[J] = ms[v] ;
/*
   fprintf(stdout, "\n vertex %d, front %d, map %d", v, J, map[J]) ;
*/
}
ndom = 0 ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) {
/*
   fprintf(stdout, "\n J = %d", J) ;
*/
   if ( map[J] != 0 ) {
      if ( (K = par[J]) != -1 ) {
         if ( map[K] == 0 ) {
            map[J] = ++ndom ;
         } else {
            map[J] = map[K] ;
         }
      } else {
         map[J] = ++ndom ;
      }
/*
      fprintf(stdout, ", in domain %d", map[J]) ;
   } else {
      fprintf(stdout, ", schur complement front") ;
*/
   }
}
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   ms[v] = map[J] ;
/*
   fprintf(stdout, "\n vertex %d, front %d, region %d", v, J, map[J]) ;
*/
}
/*
   --------------------------------------------------
   set sizes of the nvtxIV, nzfIV and opsV vectors
   to hold the domain and schur complement statistics
   --------------------------------------------------
*/
IV_setSize(nvtxIV, ndom+1) ;
IV_setSize(nzfIV,  ndom+1) ;
DV_setSize(opsDV,  ndom+1) ;
nvtxreg = IV_entries(nvtxIV) ;
nzfreg  = IV_entries(nzfIV) ;
opsreg  = DV_entries(opsDV) ;
IVzero(ndom+1, nvtxreg) ;
IVzero(ndom+1, nzfreg) ;
DVzero(ndom+1, opsreg) ;
/*
   ---------------------------------
   get the statistics for the fronts
   ---------------------------------
*/
nodwghts = ETree_nodwghts(frontETree) ;
nzfvecIV = ETree_factorEntriesIV(frontETree, symflag) ;
nzfvec   = IV_entries(nzfvecIV) ;
opsvecDV = ETree_forwardOps(frontETree, type, symflag) ;
opsvec   = DV_entries(opsvecDV) ;
/*
   ----------------------------------
   fill the region statistics vectors
   ----------------------------------
*/
for ( J = 0 ; J < nfront ; J++ ) {
   reg = map[J] ;
   nvtxreg[reg] += nodwghts[J] ;
   nzfreg[reg]  += nzfvec[J] ;
   opsreg[reg]  += opsvec[J] ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IV_free(nzfvecIV) ;
DV_free(opsvecDV) ;
IVfree(map) ;

return ; }

/*--------------------------------------------------------------------*/
