/*  misc.C  */

#include "../EGraph.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   make an element graph for a n1 x n2 grid with ncomp components

   created -- 95nov03, cca
   --------------------------------------------------------------
*/
EGraph *
EGraph_make9P ( 
   int   n1, 
   int   n2, 
   int   ncomp 
) {
EGraph    *egraph ;
int       eid, icomp, ij, ielem, jelem, m, nelem, nvtx ;
int       *list ;
#if MYDEBUG > 0
fprintf(stdout, "\n inside EGraph_make9P(%d,%d,%d)", n1, n2, ncomp) ;
fflush(stdout) ;
#endif
/*
   ---------------
   check the input
   ---------------
*/
if ( n1 <= 0 || n2 <= 0 || ncomp <= 0 ) {
   fprintf(stderr, "\n fatal error in EGraph_make9P(%d,%d,%d)"
           "\n bad input\n", n1, n2, ncomp) ;
   exit(-1) ;
}
/*
   -----------------
   create the object
   -----------------
*/
nelem = (n1 - 1)*(n2 - 1) ;
nvtx  = n1*n2*ncomp ;
egraph = EGraph_new() ;
if ( ncomp == 1 ) {
   EGraph_init(egraph, 0, nelem, nvtx, IVL_CHUNKED) ;
} else {
   EGraph_init(egraph, 1, nelem, nvtx, IVL_CHUNKED) ;
   IVfill(nvtx, egraph->vwghts, ncomp) ;
}
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
list = IVinit(4*ncomp, -1) ;
for ( jelem = 0 ; jelem < n2 - 1 ; jelem++ ) {
   for ( ielem = 0 ; ielem < n1 - 1 ; ielem++ ) {
      eid   = ielem + jelem * (n1 - 1) ;
      m  = 0 ;
      ij = ncomp*(ielem + jelem*n1) ;
      for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
         list[m++] = ij++ ;
      }
      ij = ncomp*(ielem + 1 + jelem*n1) ;
      for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
         list[m++] = ij++ ;
      }
      ij = ncomp*(ielem + (jelem+1)*n1) ;
      for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
         list[m++] = ij++ ;
      }
      ij = ncomp*(ielem + 1 + (jelem+1)*n1) ;
      for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
         list[m++] = ij++ ;
      }
      IVqsortUp(m, list) ;
      IVL_setList(egraph->adjIVL, eid, m, list) ;
   }
}
IVfree(list) ;

return(egraph) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   make an element graph for a n1 x n2 x n3 grid with ncomp components

   created -- 95nov03, cca
   -------------------------------------------------------------------
*/
EGraph *
EGraph_make27P ( 
   int   n1, 
   int   n2, 
   int   n3, 
   int   ncomp 
) {
EGraph    *egraph ;
int       eid, icomp, ijk, ielem, jelem, kelem, m, nelem, nvtx ;
int       *list ;
/*
   ---------------
   check the input
   ---------------
*/
if ( n1 <= 0 || n2 <= 0 || n3 <= 0 || ncomp <= 0 ) {
   fprintf(stderr, "\n fatal error in EGraph_make27P(%d,%d,%d,%d)"
           "\n bad input\n", n1, n2, n3, ncomp) ;
   exit(-1) ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n inside EGraph_make27P(%d,%d,%d,%d)", 
        n1, n2, n3, ncomp) ;
fflush(stdout) ;
#endif
/*
   -----------------
   create the object
   -----------------
*/
nelem = (n1 - 1)*(n2 - 1)*(n3 - 1) ;
nvtx  = n1*n2*n3*ncomp ;
egraph = EGraph_new() ;
if ( ncomp == 1 ) {
   EGraph_init(egraph, 0, nelem, nvtx, IVL_CHUNKED) ;
} else {
   EGraph_init(egraph, 1, nelem, nvtx, IVL_CHUNKED) ;
   IVfill(nvtx, egraph->vwghts, ncomp) ;
}
/*
   ----------------------------
   fill the adjacency structure
   ----------------------------
*/
list = IVinit(8*ncomp, -1) ;
for ( kelem = 0 ; kelem < n3 - 1 ; kelem++ ) {
   for ( jelem = 0 ; jelem < n2 - 1 ; jelem++ ) {
      for ( ielem = 0 ; ielem < n1 - 1 ; ielem++ ) {
         eid   = ielem + jelem*(n1-1) + kelem*(n1-1)*(n2-1);
         m   = 0 ;
         ijk = ncomp*(ielem + jelem*n1 + kelem*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + 1 + jelem*n1 + kelem*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + (jelem+1)*n1 + kelem*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + 1 + (jelem+1)*n1 + kelem*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + jelem*n1 + (kelem+1)*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + 1 + jelem*n1 + (kelem+1)*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + (jelem+1)*n1 + (kelem+1)*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         ijk = ncomp*(ielem + 1 + (jelem+1)*n1 + (kelem+1)*n1*n2) ;
         for ( icomp = 0 ; icomp < ncomp ; icomp++ ) {
            list[m++] = ijk++ ;
         }
         IVqsortUp(m, list) ;
         IVL_setList(egraph->adjIVL, eid, m, list) ;
      }
   }
}
IVfree(list) ;

return(egraph) ; }

/*--------------------------------------------------------------------*/
