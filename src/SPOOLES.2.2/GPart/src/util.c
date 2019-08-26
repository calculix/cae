/*  util.c  */

#include "../GPart.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   set the component weights from the compids[] vector

   created  -- 95oct05, cca
   modified -- 95nov29, cca
   ----------------------------------------------------
*/
void
GPart_setCweights (
   GPart   *gpart
) {
Graph   *g ;
int     ierr, ii, last, ncomp, now, nvtx, u, usize, v, w ;
int     *compids, *cweights, *list, *uadj, *vwghts ;
/*
   --------------
   check the data
   --------------
*/
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_setCweights(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
if ( (nvtx = gpart->nvtx) <= 0 || (g = gpart->g) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_setCweights(%p)"
           "\n bad Gpart object\n", gpart) ;
   exit(-1) ;
}
/*
   ----------------------------------------------------------
   set the component id of all non-multisector vertices to -1
   ----------------------------------------------------------
*/
compids = IV_entries(&gpart->compidsIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] != 0 ) {
      compids[v] = -1 ;
   }
}
/*
   ----------------------------------------------------------
   compute the number of components and set the component ids
   ----------------------------------------------------------
*/
list = IVinit(nvtx, -1) ;
ncomp = 0 ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == -1 ) {
      compids[v] = ++ncomp ;
      now = last = 0 ;
      list[now] = v ;
      while ( now <= last ) {
         u = list[now++] ;
         Graph_adjAndSize(g, u, &usize, &uadj) ;
         for ( ii = 0 ; ii < usize ; ii++ ) {
            if ( (w = uadj[ii]) < nvtx && compids[w] == -1 ) {
               compids[w] = ncomp ;
               list[++last] = w ;
            }
         }
      }
   }
}
/*
   ----------------------------
   set the number of components
   ----------------------------
*/
gpart->ncomp = ncomp ;
/*
   -------------------------
   set the component weights
   -------------------------
*/
IV_setSize(&gpart->cweightsIV, 1 + ncomp) ;
cweights = IV_entries(&gpart->cweightsIV) ;
IVzero(1 + ncomp, cweights) ;
if ( (vwghts = gpart->g->vwghts) != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      cweights[compids[v]] += vwghts[v] ;
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      cweights[compids[v]]++ ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the number of bytes taken by the object

   created -- 95oct05, cca
   ----------------------------------------------
*/
int 
GPart_sizeOf (
   GPart   *gpart
) {
int   nbytes ;

if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_sizeOf(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
nbytes = sizeof(struct _GPart) ;
nbytes += IV_size(&gpart->compidsIV)  ;
nbytes += IV_size(&gpart->cweightsIV) ;
nbytes += IV_size(&gpart->vtxMapIV)   ;

return(nbytes) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   return 1 if vertex is adjacent to only one domain
            and fill *pdomid with the domain's id
   return 0 otherwise
 
   created -- 95oct19, cca
   -------------------------------------------------
*/
int
GPart_vtxIsAdjToOneDomain (
   GPart   *gpart,
   int     v,
   int     *pdomid
) {
Graph   *g ;
int     domid, ii, nvtx, u, Vi, vsize ;
int     *compids, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL || v < 0 || (nvtx = gpart->nvtx) <= v 
   || (g = gpart->g) == NULL || pdomid == NULL ) {
   fprintf(stderr, 
           "\n fatal error in GPart_vtxIsAdjToOneDomain(%p,%d,%p)"
           "\n bad input\n", gpart, v, pdomid) ;
   exit(-1) ;
}
compids = IV_entries(&gpart->compidsIV) ;
/*
   ------------------------------------------
   fill domids[] with ids of adjacent domains
   ------------------------------------------
*/ 
Graph_adjAndSize(g, v, &vsize, &vadj) ;
domid = *pdomid = -1 ;
for ( ii = 0 ; ii < vsize ; ii++ ) {
   if ( (u = vadj[ii]) < nvtx && (Vi = compids[u]) > 0 ) {
      if ( domid == -1 ) {
         *pdomid = domid = Vi ;
      } else if ( Vi != domid ) {
         return(0) ;
      }
   }
}
if ( domid == -1 ) {
   return(0) ;
} else {
   return(1) ; }
}

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   return 1 if the partition has a valid vertex separator
          0 otherwise

   created -- 95oct18, cca
   ------------------------------------------------------
*/
int
GPart_validVtxSep (
   GPart   *gpart
) {
Graph   *g ;
int     icomp, ii, nvtx, v, vsize, w ;
int     *compids, *vadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_validVtxSep(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
nvtx    = gpart->nvtx ;
g       = gpart->g    ;
compids = IV_entries(&gpart->compidsIV) ;
/*
   ---------------------------------------------------
   loop over the vertices
   check that each non-separator vertex is adjacent to
   vertices only in its component or in the separator
   ---------------------------------------------------
*/
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( (icomp = compids[v]) != 0 ) {
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         if ( (w = vadj[ii]) < nvtx ) {
            if ( compids[w] != 0 && compids[w] != icomp ) {
               fprintf(stderr, 
"\n vertex %d, component %d, is adjacent to vertex %d, component %d",
               v, icomp, w, compids[w]) ;
               return(0) ;
            }
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   return an IV object filled with the
   weights of the component's boundaries

   created -- 96oct21, cca
   -------------------------------------
*/
IV *
GPart_bndWeightsIV (
   GPart   *gpart 
) {
Graph   *graph ;
int     icomp, ii, ncomp, nvtx, v, vsize, vwght, w ;
int     *bnd, *compids, *cweights, *mark, *vadj, *vwghts ;
IV      *bndIV ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (graph = gpart->g) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_bndWeightsIV(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
nvtx     = gpart->nvtx  ;
ncomp    = gpart->ncomp ;
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
vwghts   = graph->vwghts ;
bndIV    = IV_new() ;
IV_init(bndIV, 1 + ncomp, NULL) ;
IV_fill(bndIV, 0) ;
bnd  = IV_entries(bndIV) ;
mark = IVinit(ncomp+1, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == 0 ) {
      vwght = (vwghts == NULL) ? 1 : vwghts[v] ;
      Graph_adjAndSize(graph, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( (icomp = compids[w]) != 0 && mark[icomp] != v ) {
            mark[icomp] = v ;
            bnd[icomp] += vwght ;
         }
      }
   }
}
IVfree(mark) ;

return(bndIV) ; }

/*--------------------------------------------------------------------*/
