/*  cleanReachSet.c  */

#include "../MSMD.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   clean the vertices in the reach set

   for each v in reach set
      clean subtree list
      clean edge list
   end for

   created -- 95nov08, cca
   -------------------------------------------------
*/
void
MSMD_cleanReachSet ( 
   MSMD       *msmd,
   MSMDinfo   *info
) {
int       k, nreach ;
int       *reach ;
MSMDvtx   *v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || info == NULL ) {
   fprintf(stderr, "\n inside MSMD_cleanReachSet(%p,%p)"
           "\n bad input\n", msmd, info) ;
   exit(-1) ;
}
nreach = IV_size(&msmd->reachIV) ;
reach  = IV_entries(&msmd->reachIV) ;
/*
nreach = msmd->nreach ;
reach  = msmd->reach  ;
*/
if (  nreach < 0 || nreach > msmd->nvtx || reach == NULL ) {
   fprintf(stderr, "\n inside MSMD_cleanReachSet(%p)"
           "\n nreach = %d, reach = %p, nvtx = %d\n", 
           msmd, nreach, reach, msmd->nvtx) ;
   exit(-1) ;
}
if ( info->msglvl >= 5 ) {
   fprintf(info->msgFile, "\n inside MSMD_cleanReachSet(%p)", msmd) ;
   fflush(info->msgFile) ;
}
/*
   ---------------------------------------------------------
   clean the subtree lists of the vertices in the reach list
   ---------------------------------------------------------
*/
#if MYDEBUG > 0
   fprintf(stdout, "\n nreach = %d", nreach) ;
   for ( k = 0 ; k < nreach ; k++ ) {
      v = msmd->vertices + reach[k] ;
      fprintf(stdout, "\n <%d, %c, %c>", v->id, v->status, v->mark) ;
   }
   fflush(stdout) ;
#endif
for ( k = 0 ; k < nreach ; k++ ) {
   v = msmd->vertices + reach[k] ;
#if MYDEBUG > 1
   fprintf(stdout, "\n calling MSMD_cleanSubtreeList(%p,%d)", 
           msmd, v->id) ;
   fflush(stdout) ;
#endif
   MSMD_cleanSubtreeList(msmd, v, info) ;
}
/*
   ------------------------------------------------------
   clean the edge lists of the vertices in the reach list
   ------------------------------------------------------
*/
for ( k = 0 ; k < nreach ; k++ ) {
   v = msmd->vertices + reach[k] ;
#if MYDEBUG > 1
   fprintf(stdout, "\n calling MSMD_cleanEdgeList(%p,%d)", 
           msmd, v->id) ;
   fflush(stdout) ;
#endif
   MSMD_cleanEdgeList(msmd, v, info) ;
}

if ( info->msglvl > 3 ) {
   for ( k = 0 ; k < nreach ; k++ ) {
      v = msmd->vertices + reach[k] ;
      MSMDvtx_print(v, info->msgFile) ;
   }
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------
   clean v's subtree list of children

   created -- 95nov08, cca
   ----------------------------------
*/
void
MSMD_cleanSubtreeList ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) {
int       uid ;
IP        *ip, *nextip, *prev ;
MSMDvtx   *u ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n inside MSMD_cleanSubtreeList(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
if ( info->msglvl > 4 && info->msgFile != NULL ) {
   fprintf(info->msgFile, 
           "\n inside MSMD_cleanSubtreeList(%d)", v->id) ;
   fflush(info->msgFile) ;
}
/*
   ----------------------
   clean the subtree list
   ----------------------
*/
ip = v->subtrees ;
v->subtrees = prev = NULL ;
while ( ip != NULL ) {
   nextip = ip->next ;
   uid    = ip->val  ;
   u      = msmd->vertices + uid ;
   if ( u->par == NULL ) {
/*
      ----------------------------------------------
      adjacent subtree is still a root, keep on list
      ----------------------------------------------
*/
      if ( prev == NULL ) {
         v->subtrees = ip ;
      } else {
         prev->next = ip ;
      }
      prev = ip ;
      ip->next = NULL ;
   } else {
/*
      ----------------------------------------------------
      adjacent subtree is no longer a root, drop from list
      ----------------------------------------------------
*/
      ip->val     = -1          ;
      ip->next    = msmd->freeIP ;
      msmd->freeIP = ip          ;
   }
   ip = nextip ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   for each uncovered (v,w)
      if v->subtrees \cap w->subtrees != emptyset
         remove (v,w) from uncovered edges
      end if
   end for

   created -- 95nov08, cca
   ----------------------------------------------
*/
void
MSMD_cleanEdgeList ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) {
int       i, ierr, j, nedge, wid ;
int       *edges ;
IP        *ip1, *ip2 ;
MSMDvtx   *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n inside MSMD_cleanEdgeList(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
/*
   --------------------------------------------
   remove covered edges from the uncovered list
   --------------------------------------------
*/
nedge = v->nadj ;
edges = v->adj  ;
if ( info->msglvl > 5 ) {
   fprintf(info->msgFile, "\n inside MSMD_cleanEdgeList(%p,%p)"
           "\n %d's edges :", msmd, v, v->id) ;
   IVfp80(info->msgFile, nedge, edges, 12, &ierr) ;
   fflush(info->msgFile) ;
}
i = 0 ; j = nedge - 1 ;
while ( i <= j ) {
   wid = edges[i] ;
   w = msmd->vertices + wid ;
   if ( info->msglvl > 5 ) {
      fprintf(info->msgFile, "\n   <%d,%c>", wid, w->status) ;
      fflush(info->msgFile) ;
   }
   if ( w == v ) {
/*
      --------------------------------
      purge v from uncovered edge list
      --------------------------------
*/
      edges[i] = edges[j] ;
      edges[j] = wid ;
      j-- ;
   } else {
      switch ( w->status ) {
      case 'I' :
      case 'L' :
      case 'E' :
/*
         --------------------------------
         purge w from uncovered edge list
         --------------------------------
*/
         edges[i] = edges[j] ;
         edges[j] = wid ;
         j-- ;
         break ;
      default :
         ip1 = v->subtrees ;
         ip2 = w->subtrees ;
         if ( info->msglvl > 5 ) {
            fprintf(info->msgFile, "\n subtree list for %d :", v->id) ;
            IP_fp80(info->msgFile, ip1, 30) ;
            fprintf(info->msgFile, "\n subtree list for %d :", w->id) ;
            IP_fp80(info->msgFile, ip2, 30) ;
         }
         while ( ip1 != NULL && ip2 != NULL ) {
            if ( ip1->val > ip2->val ) {
               ip1 = ip1->next ;
            } else if ( ip1->val < ip2->val ) {
               ip2 = ip2->next ;
            } else {
/*
               ---------------------------------
               this edge has been covered, break
               ---------------------------------
*/
               edges[i] = edges[j] ;
               edges[j] = wid ;
               j-- ;
               break ;
            }
         }
         if ( ip1 == NULL || ip2 == NULL ) {
/*
            -------------------------------
            this edge was not covered, keep
            -------------------------------
*/
            i++ ;
         }
      }
   }
}
v->nadj = j + 1 ;
if ( info->msglvl > 5 ) {
   fprintf(info->msgFile, "\n leaving MSMD_cleanEdgeList(%p,%p)"
           "\n %d's edges :", msmd, v, v->id) ;
   IVfp80(info->msgFile, v->nadj, edges, 12, &ierr) ;
   fflush(info->msgFile) ;
}

return ; }

/*--------------------------------------------------------------------*/
