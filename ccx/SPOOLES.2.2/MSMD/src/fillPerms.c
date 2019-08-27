/*  fillPerms.c  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------
   fill the two permutation vectors

   created -- 96feb24, cca
   --------------------------------
*/
void
MSMD_fillPerms (
   MSMD   *msmd,
   IV     *newToOldIV,
   IV     *oldToNewIV
) {
int       count, front, iv, nfront, nvtx, pfront, root, vfront, wfront ;
int       *fch, *head, *link, *newToOld, *oldToNew, *par, *sib, 
          *vToFront ;
MSMDvtx   *v, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || (oldToNewIV == NULL && newToOldIV == NULL) ) {
   fprintf(stderr, "\n fatal error in MSMD_fillPerms(%p,%p,%p)"
           "\n bad input\n", msmd, newToOldIV, oldToNewIV) ;
   exit(-1) ;
}
nvtx = msmd->nvtx ;
if ( newToOldIV != NULL ) {
   if ( IV_size(newToOldIV) < nvtx ) {
      IV_setSize(newToOldIV, nvtx) ;  
   }
   newToOld = IV_entries(newToOldIV) ;
} else {
   newToOld = NULL ;
}
if ( oldToNewIV != NULL ) {
   if ( IV_size(oldToNewIV) < nvtx ) {
      IV_setSize(oldToNewIV, nvtx) ;  
   }
   oldToNew = IV_entries(oldToNewIV) ;
} else {
   oldToNew = NULL ;
}
/*
   -------------------------------------------
   compute the map from the vertices to fronts
   -------------------------------------------
*/
nfront = 0 ;
vToFront = IVinit(nvtx, -1) ;
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
   switch ( v->status ) {
   case 'L' :
   case 'E' :
      vToFront[iv] = nfront++ ;
      break ;
   default :
      break ;
   }
}
/*
   ------------------------
   allocate working storage
   ------------------------
*/
par  = IVinit(nfront, -1) ;
fch  = IVinit(nfront, -1) ;
sib  = IVinit(nfront, -1) ;
head = IVinit(nfront, -1) ;
link = IVinit(nvtx,   -1) ;
root = -1 ;
for ( iv = 0, v = msmd->vertices ; iv < nvtx ; iv++, v++ ) {
   switch ( v->status ) {
   case 'I' :
      w = v ;
      while ( w->status == 'I' ) {
         w = w->par ;
      }
      wfront = vToFront[w->id] ;
      link[iv] = head[wfront] ;
      head[wfront] = iv ;
      break ;
   case 'E' :
   case 'L' :
      vfront = vToFront[iv] ;
      link[iv] = head[vfront] ;
      head[vfront] = iv ;
      if ( v->par != NULL ) {
         pfront      = vToFront[v->par->id] ;
         par[vfront] = pfront ;
         sib[vfront] = fch[pfront] ;
         fch[pfront] = vfront ;
      } else {
         sib[vfront] = root ;
         root        = vfront ;
      }
      break ;
   default :
      fprintf(stderr, "\n fatal error in MSMD_fillPerms(%p,%p,%p)"
              "\n v = %d, status = %c",
              msmd, oldToNew, newToOld, v->id, v->status) ;
      fprintf(stderr, "\n vertex %d, status %c", v->id, v->status) ;
      exit(-1) ;
   }
}
#if MYDEBUG > 0
{ int   ierr ;

fprintf(stdout, "\n %d fronts, root = %d", nfront, root) ;
fprintf(stdout, "\n front par[]") ;
IVfp80(stdout, nfront, par, 80, &ierr) ;
fprintf(stdout, "\n front fch[]") ;
IVfp80(stdout, nfront, fch, 80, &ierr) ;
fprintf(stdout, "\n front sib[]") ;
IVfp80(stdout, nfront, sib, 80, &ierr) ;
fprintf(stdout, "\n head[]") ;
IVfp80(stdout, nfront, head, 80, &ierr) ;
fprintf(stdout, "\n link[]") ;
IVfp80(stdout, nvtx, head, 80, &ierr) ;
fflush(stdout) ;
}
#endif
/*
   ----------------------------------------------------------
   do a post-order traversal and fill the permutation vectors
   ----------------------------------------------------------
*/
count = 0 ;
front = root ;
while ( front != -1 ) {
   while ( fch[front] != -1 ) {
      front = fch[front] ;
   }
/*
   ---------------------------
   leaf front, number vertices
   ---------------------------
*/
   for ( iv = head[front] ; iv != -1 ; iv = link[iv] ) {
      if ( newToOld != NULL ) {
         newToOld[count] = iv ;
      }
      if ( oldToNew != NULL ) {
         oldToNew[iv] = count++ ;
      }
   }
/*
   ---------------
   find next front
   ---------------
*/
   while ( sib[front] == -1 && par[front] != -1 ) {
      front = par[front] ;
/*
      -------------------------------
      internal front, number vertices
      -------------------------------
*/
      for ( iv = head[front] ; iv != -1 ; iv = link[iv] ) {
         if ( newToOld != NULL ) {
            newToOld[count] = iv ;
         }
         if ( oldToNew != NULL ) {
            oldToNew[iv] = count++ ;
         }
      }
   }
   front = sib[front] ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(par)      ;
IVfree(fch)      ;
IVfree(sib)      ;
IVfree(head)     ;
IVfree(link)     ;
IVfree(vToFront) ;

return ; }

/*--------------------------------------------------------------------*/
