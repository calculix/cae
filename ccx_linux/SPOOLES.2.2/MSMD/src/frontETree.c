/*  frontETree.c  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   create and return an ETree object that holds the front tree.

   created  -- 96jun23, cca
   ------------------------------------------------------------
*/
ETree *
MSMD_frontETree (
   MSMD   *msmd
) {
ETree     *etree ;
int       front, iv, nfront, nvtx, root ;
int       *bndwghts, *fch, *nodwghts, *par, *sib, *vtxToFront ;
MSMDvtx   *v, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_frontETree(%p)"
           "\n bad input\n", msmd) ;
   exit(-1) ;
}
nvtx = msmd->nvtx ;
/*
   --------------------------
   count the number of fronts
   --------------------------
*/
nfront = 0 ;
fch = IVinit(nvtx, -1) ;
sib = IVinit(nvtx, -1) ;
root = -1 ;
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n vertex %d, status %c, wght %d",
           v->id, v->status, v->wght) ;
/*
   MSMDvtx_print(v, stdout) ;
*/
   fflush(stdout) ;
#endif
   switch ( v->status ) {
   case 'L' :
   case 'E' :
      if ( (w = v->par) != NULL ) {
         sib[v->id] = fch[w->id] ;
         fch[w->id] = v->id ;
      } else {
         sib[v->id] = root ;
         root = v->id ;
      }
#if MYDEBUG > 0
   fprintf(stdout, ", new front %d", nfront) ;
   fflush(stdout) ;
#endif
      nfront++ ;
      break ;
   default :
      break ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n %d fronts", nfront) ;
fflush(stdout) ;
#endif
/*
   ---------------------------
   initialize the ETree object
   ---------------------------
*/
etree = ETree_new() ;
ETree_init1(etree, nfront, nvtx) ;
nodwghts   = IV_entries(etree->nodwghtsIV) ;
bndwghts   = IV_entries(etree->bndwghtsIV) ;
vtxToFront = IV_entries(etree->vtxToFrontIV) ;
/*
   ----------------------------------------------
   fill the vtxToFront[] vector so representative 
   vertices are mapped in a post-order traversal
   ----------------------------------------------
*/
nfront = 0 ;
iv = root ;
while ( iv != -1 ) {
   while ( fch[iv] != -1 ) {
       iv = fch[iv] ;
   }
   v = msmd->vertices + iv ;
   vtxToFront[iv] = nfront++ ;
#if MYDEBUG > 0
   fprintf(stdout, "\n v = %d, vwght = %d, vtxToFront[%d] = %d", 
           v->id, v->wght, iv, vtxToFront[iv]) ;
   fflush(stdout) ;
#endif
   while ( sib[iv] == -1 && v->par != NULL ) {
      v = v->par ;
      iv = v->id ;
      vtxToFront[iv] = nfront++ ;
#if MYDEBUG > 0
      fprintf(stdout, "\n v = %d, vwght = %d, vtxToFront[%d] = %d", 
              v->id, v->wght, iv, vtxToFront[iv]) ;
      fflush(stdout) ;
#endif
   }
   iv = sib[iv] ;
}
IVfree(fch) ;
IVfree(sib) ;
/*
   --------------------------------------------------------------
   fill in the vertex-to-front map for indistinguishable vertices
   --------------------------------------------------------------
*/
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n v %d, wght = %d, status %c", 
           v->id, v->wght, v->status) ;
   fflush(stdout) ;
#endif
   switch ( v->status ) {
   case 'I' :
#if MYDEBUG > 0
      fprintf(stdout, "\n I : v %d", v->id) ;
      fflush(stdout) ;
#endif
      w = v ;
      while ( w->par != NULL && w->status == 'I' ) {
         w = w->par ;
#if MYDEBUG > 0
/*
         fprintf(stdout, " --> %d", w->id) ;
*/
         fprintf(stdout, " %d", w->id) ;
         fflush(stdout) ;
#endif
      }
#if MYDEBUG > 0
      fprintf(stdout, ", w %d, status %c", w->id, w->status) ;
      fflush(stdout) ;
#endif
      switch ( w->status ) {
      case 'L' :
      case 'E' :
         vtxToFront[v->id] = vtxToFront[w->id] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n I: vtxToFront[%d] = %d", 
                 iv, vtxToFront[iv]) ;
         fflush(stdout) ;
#endif
         break ;
      default :
#if MYDEBUG > 0
         fprintf(stdout, "\n wow, v->rootpar = %d, status %c",
                 w->id, w->status) ;
         fflush(stdout) ;
#endif
         break ;
      }
   }
}
/*
   ------------------------------------------------------------
   now fill in the parent Tree field, node and boundary weights 
   ------------------------------------------------------------
*/
par = etree->tree->par ;
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n v %d, status %c", v->id, v->status) ;
   fflush(stdout) ;
#endif
   switch ( v->status ) {
   case 'L' :
   case 'E' :
      front = vtxToFront[iv] ;
#if MYDEBUG > 0
      fprintf(stdout, ", front %d", front) ;
      fflush(stdout) ;
#endif
      if ( (w = v->par) != NULL ) {
         par[vtxToFront[v->id]] = vtxToFront[w->id] ;
#if MYDEBUG > 0
         fprintf(stdout, ", par[%d] = %d", front, par[front]) ;
         fflush(stdout) ;
#endif
      }
      bndwghts[front] = v->bndwght ;
      nodwghts[front] = v->wght    ;
      break ;
   default :
      break ;
   }
}
/*
   -------------------------
   set the other tree fields
   -------------------------
*/
Tree_setFchSibRoot(etree->tree) ;

return(etree) ; }

/*--------------------------------------------------------------------*/
