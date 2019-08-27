/*  makeSchurComplement.c  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose -- 

   if the elimination has halted before all the stages have been 
   eliminated, then create the schur complement graph and the map 
   from the original vertices those in the schur complement graph.

   schurGraph -- Graph object to contain the schur complement graph
   VtoPhi     -- IV object to contain the map from vertices in V
                 to schur complement vertices in Phi

   created -- 97feb01, cca
   ----------------------------------------------------------------
*/
void
MSMD_makeSchurComplement (
   MSMD    *msmd,
   Graph   *schurGraph,
   IV      *VtoPhiIV
) {
int       nedge, nPhi, nvtx, totewght, totvwght ;
int       *mark, *rep, *VtoPhi, *vwghts ;
int       count, *list ;
int       ierr, ii, size, *adj ;
int       phi, psi, tag ;
IP        *ip ;
IVL       *adjIVL ;
MSMDvtx   *u, *v, *vertices, *vfirst, *vlast, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || schurGraph == NULL || VtoPhiIV == NULL ) {
   fprintf(stderr, 
           "\n\n fatal error in MSMD_makeSchurComplement(%p,%p,%p)"
           "\n bad input\n", msmd, schurGraph, VtoPhiIV) ;
   exit(-1) ;
}
vertices = msmd->vertices ;
nvtx     = msmd->nvtx     ;
/*
   -------------------------------------
   initialize the V-to-Phi map IV object
   -------------------------------------
*/
IV_clearData(VtoPhiIV) ;
IV_setSize(VtoPhiIV, nvtx) ;
IV_fill(VtoPhiIV, -2) ;
VtoPhi = IV_entries(VtoPhiIV) ;
/*
   ---------------------------------------------
   count the number of Schur complement vertices
   ---------------------------------------------
*/
vfirst = vertices ;
vlast  = vfirst + nvtx - 1 ;
nPhi   = 0 ;
for ( v = vfirst ; v <= vlast ; v++ ) {
#if MYDEBUG > 0
   fprintf(stdout, "\n v->id = %d, v->status = %c", v->id, v->status) ;
   fflush(stdout) ;
#endif
   switch ( v->status ) {
   case 'L' :
   case 'E' :
   case 'I' :
      break ;
   case 'B' :
      VtoPhi[v->id] = nPhi++ ;
#if MYDEBUG > 0
      fprintf(stdout, ", VtoPhi[%d] = %d", v->id, VtoPhi[v->id]) ;
      fflush(stdout) ;
#endif
      break ;
   default :
      break ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n nPhi = %d", nPhi) ;
fflush(stdout) ;
#endif
/*
   ----------------------------------------------------
   get the representative vertex id for each Phi vertex
   ----------------------------------------------------
*/
rep = IVinit(nPhi, -1) ;
for ( v = vfirst ; v <= vlast ; v++ ) {
   if ( (phi = VtoPhi[v->id]) >= 0 ) {
#if MYDEBUG > 0
      fprintf(stdout, "\n rep[%d] = %d", phi, v->id) ;
      fflush(stdout) ;
#endif
      rep[phi] = v->id ;
   }
}
/*
   ------------------------------------------
   set the map for indistinguishable vertices
   ------------------------------------------
*/
for ( v = vfirst ; v <= vlast ; v++ ) {
   if ( v->status == 'I' ) {
      w = v ;
      while ( w->status == 'I' ) {
         w = w->par ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n v = %d, status = %c, w = %d, status = %c", 
              v->id, v->status, w->id, w->status) ;
      fflush(stdout) ;
#endif
      VtoPhi[v->id] = VtoPhi[w->id] ;
   }
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n VtoPhi") ;
IV_writeForHumanEye(VtoPhiIV, stdout) ;
fflush(stdout) ;
#endif
/*
   ---------------------------
   initialize the Graph object
   ---------------------------
*/
Graph_clearData(schurGraph) ;
Graph_init1(schurGraph, 1, nPhi, 0, 0, IVL_CHUNKED, IVL_CHUNKED) ;
adjIVL = schurGraph->adjIVL ;
vwghts = schurGraph->vwghts ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n schurGraph initialized, nvtx = %d",
        schurGraph->nvtx) ;
fflush(stdout) ;
#endif
/*
   -------------------------------
   fill the vertex adjacency lists
   -------------------------------
*/
mark = IVinit(nPhi, -1) ;
list = IVinit(nPhi, -1) ;
nedge = totvwght = totewght = 0 ;
for ( phi = 0 ; phi < nPhi ; phi++ ) {
/*
   -----------------------------
   get the representative vertex
   -----------------------------
*/
   v = vfirst + rep[phi] ; 
#if MYDEBUG > 0
   fprintf(stdout, "\n phi = %d, v = %d", phi, v->id) ;
   fflush(stdout) ;
   MSMDvtx_print(v, stdout) ;
   fflush(stdout) ;
#endif
   count = 0 ;
   tag   = v->id ;
/*
   ---------------------------
   load self in adjacency list
   ---------------------------
*/
   mark[phi] = tag ;
   totewght += v->wght * v->wght ;
#if MYDEBUG > 0
   fprintf(stdout, "\n    mark[%d] = %d", phi, mark[phi]) ;
   fflush(stdout) ;
#endif
   list[count++] = phi ;
/*
   ----------------------------------------
   load boundary lists of adjacent subtrees 
   ----------------------------------------
*/
   for ( ip = v->subtrees ; ip != NULL ; ip = ip->next ) {
      u    = vertices + ip->val ;
      size = u->nadj ;
      adj  = u->adj  ;
#if MYDEBUG > 0
      fprintf(stdout, "\n    subtree %d :", u->id) ;
      IVfp80(stdout, size, adj, 15, &ierr) ;
      fflush(stdout) ;
#endif
      for ( ii = 0 ; ii < size ; ii++ ) {
         w = vertices + adj[ii] ;
#if MYDEBUG > 0
         fprintf(stdout, "\n       w %d, status %c, psi %d",
                 w->id, w->status, VtoPhi[w->id]) ;
         fflush(stdout) ;
#endif
         if ( (psi = VtoPhi[w->id]) != -2 && mark[psi] != tag ) {
            mark[psi] = tag ;
#if MYDEBUG > 0
            fprintf(stdout, ", mark[%d] = %d", psi, mark[psi]) ;
            fflush(stdout) ;
#endif
            list[count++] = psi ;
            totewght += v->wght * w->wght ;
         }
      }
   }
/*
   ----------------------
   load adjacent vertices 
   ----------------------
*/
   size = v->nadj ;
   adj  = v->adj  ;
   for ( ii = 0 ; ii < size ; ii++ ) {
      w = vertices + adj[ii] ;
      if ( (psi = VtoPhi[w->id]) != -2 && mark[psi] != tag ) {
         mark[psi] = tag ;
         list[count++] = psi ;
         totewght += v->wght * w->wght ;
      }
   }
/*
   ---------------------------------------------
   sort the list and inform adjacency IVL object
   ---------------------------------------------
*/
   IVqsortUp(count, list) ;
   IVL_setList(adjIVL, phi, count, list) ;
/*
   --------------------------------------
   set the vertex weight and increment 
   the total vertex weight and edge count
   --------------------------------------
*/
   vwghts[phi] =  v->wght ;
   totvwght    += v->wght ;
   nedge       += count   ;
}
schurGraph->totvwght = totvwght ;
schurGraph->nedges   = nedge    ;
schurGraph->totewght = totewght ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;
IVfree(mark) ;
IVfree(rep)  ;

return ; }
 
/*--------------------------------------------------------------------*/
