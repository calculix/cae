/*  instance.c  */

#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   return the number of fronts

   created -- 97feb28, cca
   ---------------------------
*/
int
ETree_nfront (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_nfront(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->nfront) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   return the number of vertices

   created -- 97feb28, cca
   -----------------------------
*/
int
ETree_nvtx (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_nvtx(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->nvtx) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   return a pointer to the Tree object

   created -- 97feb28, cca
   -----------------------------------
*/
Tree *
ETree_tree (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_tree(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->tree) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------
   return the root of the tree

   created -- 97feb28, cca
   ---------------------------
*/
int
ETree_root (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_root(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->tree->root) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return a pointer to the parent vector

   created -- 97feb28, cca
   -------------------------------------
*/
int *
ETree_par (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_par(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->tree->par) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return a pointer to the first child vector

   created -- 97feb28, cca
   ------------------------------------------
*/
int *
ETree_fch (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_fch(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->tree->fch) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   return a pointer to the sibling vector

   created -- 97feb28, cca
   --------------------------------------
*/
int *
ETree_sib (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_sib(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->tree->sib) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return a pointer to the nodwghts IV object

   created -- 97feb28, cca
   ------------------------------------------
*/
IV *
ETree_nodwghtsIV (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_nodwghtsIV(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->nodwghtsIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   return a pointer to the nodwghts int vector

   created -- 97feb28, cca
   -------------------------------------------
*/
int *
ETree_nodwghts (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->nodwghtsIV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_nodwghts(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(IV_entries(etree->nodwghtsIV)) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   return a pointer to the bndwghts IV object

   created -- 97feb28, cca
   ------------------------------------------
*/
IV *
ETree_bndwghtsIV (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_bndwghtsIV(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->bndwghtsIV) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   return a pointer to the bndwghts int vector

   created -- 97feb28, cca
   -------------------------------------------
*/
int *
ETree_bndwghts (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->bndwghtsIV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_bndwghts(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(IV_entries(etree->bndwghtsIV)) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   return a pointer to the vtxToFront IV object

   created -- 97feb28, cca
   --------------------------------------------
*/
IV *
ETree_vtxToFrontIV (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_vtxToFrontIV(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(etree->vtxToFrontIV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   return a pointer to the vtxToFront int vector

   created -- 97feb28, cca
   ---------------------------------------------
*/
int *
ETree_vtxToFront (
   ETree   *etree
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || etree->vtxToFrontIV == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_vtxToFront(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
return(IV_entries(etree->vtxToFrontIV)) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- return the number of internal degrees 
              of freedom in front J

   created -- 97may23, cca
   ------------------------------------------------
*/
int
ETree_frontSize (
   ETree   *etree,
   int     J
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || J < 0 || J >= etree->nfront ) {
   fprintf(stderr, "\n fatal error in ETree_frontSize(%p,%d)"
           "\n bad input\n", etree, J) ;
   exit(-1) ;
}
return(etree->nodwghtsIV->vec[J]) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- return the number of external degrees 
              of freedom in front J

   created -- 97may23, cca
   ------------------------------------------------
*/
int
ETree_frontBoundarySize (
   ETree   *etree,
   int     J
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || J < 0 || J >= etree->nfront ) {
   fprintf(stderr, "\n fatal error in ETree_frontBoundarySize(%p,%d)"
           "\n bad input\n", etree, J) ;
   exit(-1) ;
}
return(etree->bndwghtsIV->vec[J]) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- compute the maximum number of indices and entries 
              in a front

   symflag = SPOOLES_SYMMETRIC or SPOOLES_HERMITIAN --> 
      count only column indices
      count upper entries in (1,1) block and (1,2) block
   symflag = SPOOLES_NONSYMMETRIC --> 
      count row and column indices
      count entries in (1,1), (1,2) and (2,1) blocks

   created -- 97may23, cca
   ------------------------------------------------------------
*/
void
ETree_maxNindAndNent (
   ETree   *etree,
   int     symflag,
   int     *pmaxnind,
   int     *pmaxnent
) {
int   J, maxnent, maxnind, nDJ, nent, nfront, nind, nUJ ;
int   *nodwghts, *bndwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_maxNindAndNent(%p,%d)"
           "\n bad input\n", etree, symflag) ;
   exit(-1) ;
}
nfront   = etree->nfront ;
nodwghts = ETree_nodwghts(etree) ;
bndwghts = ETree_bndwghts(etree) ;
for ( J = 0, maxnent = maxnind = 0 ; J < nfront ; J++ ) {
   nDJ = nodwghts[J] ;
   nUJ = bndwghts[J] ;
   switch ( symflag ) {
   case SPOOLES_SYMMETRIC :
   case SPOOLES_HERMITIAN :
      nind = nDJ + nUJ ;
      nent = (nDJ*(nDJ+1))/2 + nDJ*nUJ ;
      break ;
   case SPOOLES_NONSYMMETRIC :
      nind = 2*(nDJ + nUJ) ;
      nent = nDJ*(nDJ + 2*nUJ) ;
      break ;
   default :
      fprintf(stderr, "\n fatal error in ETree_maxNindAndNent(%p,%d)"
              "\n bad symflag\n", etree, symflag) ;
      exit(-1) ;
      break ;
   }
   if ( maxnind < nind ) {
      maxnind = nind ;
   }
   if ( maxnent < nent ) { 
      maxnent = nent ;
   }
}
*pmaxnind = maxnind ;
*pmaxnent = maxnent ;

return ; }

/*--------------------------------------------------------------------*/
