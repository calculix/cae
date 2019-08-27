/*  TwoSetViaBKL.c  */

#include "../GPart.h"
#include "../../BKL.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------------
   given a domain decomposition, find a bisector
   1. construct the domain/segment graph
   2. use block kernihan-lin to get an initial bisector

   alpha   -- cost function parameter for BKL
   seed    -- random number seed
   cpus    -- array to store CPU times
              cpus[0] -- time to find domain/segment map
              cpus[1] -- time to find domain/segment bipartite graph
              cpus[2] -- time to find two-set partition

   return value -- cost of the partition

   created  -- 96mar09, cca
   -----------------------------------------------------------------
*/
double
GPart_TwoSetViaBKL (
   GPart       *gpart,
   double      alpha,
   int         seed,
   double      cpus[]
) {
BKL      *bkl ;
BPG      *bpg ;
double   t1, t2 ;
FILE     *msgFile ;
float    bestcost ;
Graph    *g, *gc ;
int      c, flag, ierr, msglvl, ndom, nseg, nvtx, v ;
int      *compids, *cweights, *dscolors, *dsmap, *vwghts ;
IV       *dsmapIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL || cpus == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_DDsep(%p,%f,%d,%p)"
          "\n bad input\n", gpart, alpha, seed, cpus) ;
   exit(-1) ;
}
g        = gpart->g        ;
nvtx     = gpart->nvtx     ;
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
vwghts   = g->vwghts      ;
msglvl   = gpart->msglvl  ;
msgFile  = gpart->msgFile ;
/*
   HARDCODE THE ALPHA PARAMETER.
*/
alpha = 1.0 ;
/*
   ------------------------------
   (1) get the domain/segment map 
   (2) get the compressed graph
   (3) create the bipartite graph
   ------------------------------
*/
MARKTIME(t1) ;
dsmapIV = GPart_domSegMap(gpart, &ndom, &nseg) ;
dsmap = IV_entries(dsmapIV) ;
MARKTIME(t2) ;
cpus[0] = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %9.5f : generate domain-segment map",
           t2 - t1) ;
   fprintf(msgFile, "\n ndom = %d, nseg = %d", ndom, nseg) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------
   create the domain/segment bipartite graph
   -----------------------------------------
*/
MARKTIME(t1) ;
gc = Graph_compress(gpart->g, dsmap, 1) ;
bpg = BPG_new() ;
BPG_init(bpg, ndom, nseg, gc) ;
MARKTIME(t2) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %9.5f : create domain-segment graph",
           t2 - t1) ;
   fflush(msgFile) ;
}
cpus[1] = t2 - t1 ;
if ( msglvl > 2 ) {
   if ( bpg->graph->vwghts != NULL ) {
      fprintf(msgFile, "\n domain weights :") ;
      IVfp80(msgFile, bpg->nX, bpg->graph->vwghts, 17, &ierr) ;
      fprintf(msgFile, "\n segment weights :") ;
      IVfp80(msgFile, bpg->nY, bpg->graph->vwghts+bpg->nX, 18, &ierr) ;
      fflush(msgFile) ;
   }
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n dsmapIV ") ;
   IV_writeForHumanEye(dsmapIV, msgFile) ;
   fprintf(msgFile, "\n\n domain/segment bipartite graph ") ;
   BPG_writeForHumanEye(bpg, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   create and initialize the BKL object
   ------------------------------------
*/
MARKTIME(t1) ;
flag = 5 ;
bkl = BKL_new() ;
BKL_init(bkl, bpg, alpha) ;
BKL_setInitPart(bkl, flag, seed, NULL) ;
bestcost = BKL_evalfcn(bkl) ;
gpart->ncomp = 2 ;
MARKTIME(t2) ;
cpus[2] = t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %9.5f : initialize BKL object", t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n BKL : flag = %d, seed = %d", flag, seed) ;
   fprintf(msgFile, ", initial cost = %.2f", bestcost) ;
   fflush(msgFile) ;
   fprintf(msgFile, ", cweights = < %d %d %d >",
           bkl->cweights[0], bkl->cweights[1], bkl->cweights[2]) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n colors") ;
   IVfp80(msgFile, bkl->nreg, bkl->colors, 80, &ierr) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n BKL initial weights : ") ;
   IVfp80(msgFile, 3, bkl->cweights, 25, &ierr) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   improve the partition via fidmat
   --------------------------------
*/
MARKTIME(t1) ;
bestcost = BKL_fidmat(bkl) ;
MARKTIME(t2) ;
cpus[2] += t2 - t1 ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n CPU %9.5f : improve the partition via fidmat",
           t2 - t1) ;
   fflush(msgFile) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n BKL : %d passes", bkl->npass) ;
   fprintf(msgFile, ", %d flips", bkl->nflips) ;
   fprintf(msgFile, ", %d gainevals", bkl->ngaineval) ;
   fprintf(msgFile, ", %d improve steps", bkl->nimprove) ;
   fprintf(msgFile, ", cost = %9.2f", bestcost) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, 
        "\n BKL STATS < %9d %9d %9d > %9.2f < %4d %4d %4d %4d %4d >",
        bkl->cweights[0], bkl->cweights[1], bkl->cweights[2],
        bestcost, bkl->npass, bkl->npatch, bkl->nflips, bkl->nimprove,
        bkl->ngaineval) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n colors") ;
   IVfp80(msgFile, bkl->nreg, bkl->colors, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ----------------------------
   set compids[] and cweights[]
   ----------------------------
*/
MARKTIME(t1) ;
dscolors = bkl->colors ;
gpart->ncomp = 2 ;
IV_setSize(&gpart->cweightsIV, 3) ;
cweights = IV_entries(&gpart->cweightsIV) ;
cweights[0] = cweights[1] = cweights[2] = 0 ;
if ( vwghts == NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      compids[v] = c = dscolors[dsmap[v]] ;
      cweights[c]++ ;
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      compids[v] = c = dscolors[dsmap[v]] ;
      cweights[c] += vwghts[v] ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n BKL partition : < %d %d %d >",
           cweights[0], cweights[1], cweights[2]) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------
   free the BKL object, the BPG object
   and the domain/segment map IV object
   ------------------------------------
*/
BKL_free(bkl) ;
IV_free(dsmapIV) ;
BPG_free(bpg) ;
MARKTIME(t2) ;
cpus[2] += t2 - t1 ;

return((double) bestcost) ; }

/*--------------------------------------------------------------------*/
