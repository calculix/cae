/*  metrics.c  */

#include "../ETree.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   return an IV object with the weights 
   of the vertices in each front.

   created -- 96jun23, cca
   ------------------------------------
*/
IV *
ETree_nvtxMetric (
   ETree   *etree
) {
IV   *metricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || etree->nfront <= 0 || etree->nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nvtxMetric(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
metricIV = IV_new() ;
IV_init(metricIV, etree->nfront, NULL) ;
IVcopy(etree->nfront, IV_entries(metricIV), 
                      IV_entries(etree->nodwghtsIV)) ;

return(metricIV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   return an IV object with the number 
   of factor entries in each front.

   symflag -- symmetryflag 
      SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 96jun23, cca
   ---------------------------------------------------------------
*/
IV *
ETree_nentMetric (
   ETree   *etree,
   int     flag
) {
int    front, nfront, nb, nv ;
int    *bndwghts, *metric, *nodwghts ;
IV     *metricIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL 
   || (nfront = etree->nfront) <= 0 || etree->nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nentMetric(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
metricIV = IV_new() ;
IV_init(metricIV, nfront, NULL) ;
metric   = IV_entries(metricIV) ;
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
if ( flag == 1 ) {
   for ( front = 0 ; front < nfront ; front++ ) {
      nv = nodwghts[front] ;
      nb = bndwghts[front] ;
      metric[front] = (nv*(nv+1))/2 + nv*nb ;
   }
} else if ( flag == 2 ) {
   for ( front = 0 ; front < nfront ; front++ ) {
      nv = nodwghts[front] ;
      nb = bndwghts[front] ;
      metric[front] = nv*nv + 2*nv*nb ;
   }
}

return(metricIV) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   return a DV object with the number 
   of factor operations in each front.

   type -- type of entries,
      SPOOLES_REAL or SPOOLES_COMPLEX

   symflag -- symmetryflag,
      SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC

   created -- 96jun23, cca
   ---------------------------------------------------------------
*/
DV *
ETree_nopsMetric (
   ETree   *etree,
   int     type,
   int     symflag
) {
DV       *metricDV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  etree == NULL || etree->nfront <= 0 || etree->nvtx <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_nopsMetric(%p)"
           "\n bad input\n", etree) ;
   exit(-1) ;
}
metricDV = ETree_forwardOps(etree, type, symflag) ;

return(metricDV) ; }

/*--------------------------------------------------------------------*/
