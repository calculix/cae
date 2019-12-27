/*  eliminate.c  */

#include "../MSMD.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   eliminate all nodes in a stage

   created -- 96feb25, cca
   ---------------------------------------------------------------------
*/
void
MSMD_eliminateStage ( 
   MSMD       *msmd,
   MSMDinfo   *info
) {
int       ierr, ii, jj, iv, nelim, nreach, stage, step ;
int       *reach ;
IV        *reachIV ;
MSMDvtx   *v ;
/*
   ---------------
   check the input
   ---------------
*/
if (  msmd == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_eliminateStage(%p,%p)"
           "\n bad input\n", msmd, info) ;
   exit(-1) ;
}
stage = info->istage ;
/*
   -----------------------------------------------
   load the reach set with all nodes in this stage
   -----------------------------------------------
*/
reachIV = &msmd->reachIV ;
IV_setSize(reachIV, 0) ;
for ( iv = 0, v = msmd->vertices ; iv < msmd->nvtx ; iv++, v++ ) {
   if ( v->status != 'I' ) {
      if ( v->stage == stage ) {
         IV_push(reachIV, v->id) ;
         v->status = 'R' ;
      } else if ( v->stage > stage || v->stage < 0 ) {
         v->status = 'B' ;
      }
   }
}
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n after loading reach set") ;
   IV_fp80(reachIV, info->msgFile, 80, &ierr) ;
   fflush(info->msgFile) ;
}
if ( info->seed > 0 ) {
   IV_shuffle(reachIV, info->seed) ;
}
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n reach set at stage %d", stage) ;
   IV_fp80(reachIV, info->msgFile, 80, &ierr) ;
   fflush(info->msgFile) ;
}
/*
   ------------------------------------
   do an initial update of the vertices
   ------------------------------------
*/
MSMD_update(msmd, info) ;
if ( info->msglvl > 4 ) {
   fprintf(info->msgFile, "\n\n after initial update") ;
   fflush(info->msgFile) ;
}
IV_setSize(reachIV, 0) ;
/*
   -----------
   elimination
   -----------
*/
step = 0 ; 
while ( 1 ) {
   if ( info->msglvl > 1 ) {
      fprintf(info->msgFile, "\n\n ##### stage %d, elimination step %d",
              stage, step) ;
      fflush(info->msgFile) ;
   }
   nelim = MSMD_eliminateStep(msmd, info) ;
   if ( nelim == 0 ) {
      break ;
   }
/*
   -----------------------------------------
   for each node in the reach set, clean its 
   subtree list and list of uncovered edges
   -----------------------------------------
*/
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n calling MSMD_cleanReachSet()") ;
      fprintf(info->msgFile, "\n reach set") ;
      IV_fp80(reachIV, info->msgFile, 80, &ierr) ;
      fflush(info->msgFile) ;
   }
   MSMD_cleanReachSet(msmd, info) ;
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n return from MSMD_cleanReachSet()") ;
      fflush(info->msgFile) ;
   }
/*
   ------------------
   compress the graph
   ------------------
*/
   MSMD_findInodes(msmd, info) ;
/*
   ----------------------------------------------
   clean the reach set of indistinguishable nodes
   and any nodes not on the present stage
   ----------------------------------------------
*/
   nreach = IV_size(reachIV) ;
   reach  = IV_entries(reachIV) ;
   for ( ii = jj = 0 ; ii < nreach ; ii++ ) {
      if ( reach[ii] < 0 || reach[ii] >= msmd->nvtx ) {
         fprintf(stderr, "\n fatal error in MSMD_eliminateStage()"
                 "\n reach[%d] = %d", ii, reach[ii]) ;
         exit(-1) ;
      }
      v = msmd->vertices + reach[ii] ;
      if ( v->status == 'I' ) {
         continue ;
      } else if ( v->stage != stage ) {
         v->status = 'B' ;
      } else {
         reach[jj++] = v->id ;
      }
   }
   IV_setSize(reachIV, jj) ;
   if ( info->msglvl > 2 ) {
      fprintf(info->msgFile, 
              "\n\n after cleaning reach set, nreach = %d",
              IV_size(reachIV)) ;
      fprintf(info->msgFile, "\n reach :") ;
      IV_fp80(reachIV, info->msgFile, 8, &ierr) ;
      fflush(info->msgFile) ;
   }
/*
   ----------------------------------
   update the nodes on the reach set
   ----------------------------------
*/
   MSMD_update(msmd, info) ;
   if ( info->msglvl > 2 ) {
      fprintf(info->msgFile, "\n\n return from update") ;
      fflush(info->msgFile) ;
   }
   IV_setSize(reachIV, 0) ;
/*
   ------------------------------
   increment the elimination step
   ------------------------------
*/
   step++ ;
}
if ( info->msglvl > 2 ) {
   fprintf(info->msgFile, "\n stage %d over, %d steps",
           stage, step) ;
   fflush(info->msgFile) ;
}
/*
   --------------------------------------
   set the number of steps for this stage
   --------------------------------------
*/
info->stageInfo->nstep = step ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to eliminate an independent set of vertices

   created -- 95feb25, cca
   ------------------------------------------------------
*/
int
MSMD_eliminateStep ( 
   MSMD       *msmd,
   MSMDinfo   *info
) {
int       deg, maxdeg, mindeg, ncol, nelim, nrow, vid, weight ;
MSMDvtx   *v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_eliminate(%p,%p)"
           "\n bad input\n", msmd, info) ;
   exit(-1) ;
}
/*
   --------------------
   check for empty heap
   --------------------
*/
if ( msmd->heap->size == 0 ) {
   return(0) ;
}
if ( info->msglvl > 2 ) {
   fprintf(info->msgFile, "\n step %d", info->stageInfo->nstep) ;
   fflush(info->msgFile) ;
}
/*
   ----------------------------------------
   increase the number of elimination steps
   ----------------------------------------
*/
info->stageInfo->nstep++ ;
/*
   -----------------------------------------
   loop while nodes of minimum degree remain
   -----------------------------------------
*/
nelim = weight = 0 ;
IIheap_root(msmd->heap, &vid, &mindeg) ;
if ( info->stepType <= 1.0 ) {
   maxdeg = mindeg ;
} else {
   maxdeg = (int) (info->stepType*mindeg) ;
}
do {
   IIheap_root(msmd->heap, &vid, &deg) ;
   if ( deg > maxdeg ) {
      break ;
   }
   v = msmd->vertices + vid ;
/*
   ------------------------------------------------------------------
   increment the weight of this elimination step, remove v from the 
   degree heap and put the vertex on the stack of eliminated vertices
   ------------------------------------------------------------------
*/
   if ( info->msglvl > 1 ) {
      fprintf(info->msgFile, 
              "\n    eliminating vertex %d, weight %d, deg %d", 
              vid, v->wght, deg) ;
      fflush(info->msgFile) ;
   }
   info->stageInfo->nfront++ ;
   info->stageInfo->welim += v->wght ;
   nelim++ ;
   weight += v->wght ;
   IIheap_remove(msmd->heap, vid) ;
/*
   ---------------------------------------------------------------
   call MSMD_eliminateVtx(v) to update v's adjacency structure and
   clean the adjacency structures of its adjacent vertices
   ---------------------------------------------------------------
*/
   MSMD_eliminateVtx(msmd, v, info) ;
/*
   ------------------------------------------
   increment the storage and operation counts
   ------------------------------------------
*/
   ncol = v->wght    ;
   nrow = v->bndwght ;
   info->stageInfo->nfind += nrow + ncol ;
   info->stageInfo->nzf   += nrow*ncol + (ncol*(ncol+1))/2 ;
   info->stageInfo->ops   += ((double) ncol*nrow)*((double) nrow+ncol+1)
                          +  ((double) ncol)*((double) (ncol+1))
                               *((double) (2*ncol+1))/6 ;
/*
fprintf(stdout, "\n nfind = %d, nzf = %d, ops = %.0f",
        info->stageInfo->nfind, info->stageInfo->nzf,
        info->stageInfo->ops) ;
*/
   if ( info->stepType < 1.0 ) {
/*
      --------------------------------------------
      we are not using multiple elimination, break
      --------------------------------------------
*/
      break ;
   }
} while ( msmd->heap->size > 0 ) ;

return(nelim) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- eliminate vertex v
      1) create v's boundary list
      2) merge boundary list onto reach list
      3) for each vertex in the boundary
         3.1) add v to the subtree list

   created -- 96feb25, cca
   -----------------------------------------
*/
void 
MSMD_eliminateVtx ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info 
) {
int       i, ierr, j, nadj, nbnd, nedge, uid, wid, wght ;
int       *adj, *bnd, *edges ;
IP        *ip, *ip2, *prev ;
IV        *reachIV ;
MSMDvtx   *u, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_eliminateVtx(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
adj     = IV_entries(&msmd->ivtmpIV) ;
reachIV = &msmd->reachIV ;
/*
   -----------------------------
   create the boundary set for v
   -----------------------------
*/
v->mark = 'X' ;
if ( v->subtrees == NULL ) {
/*
   --------------------------------------------------------
   v is a leaf, look at its uncovered edge list, move v and 
   any indistinguishable vertices to the end of the list
   --------------------------------------------------------
*/
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n vertex %d is a leaf", v->id) ;
      fflush(info->msgFile) ;
   }
   v->status = 'L' ;
   nedge = v->nadj ;
   edges = v->adj  ;
   i = 0 ; j = nedge - 1 ;
   while ( i <= j ) {
      wid = edges[i] ;
      w   = msmd->vertices + wid ;
      if ( w == v || w->status == 'I' ) {
         edges[i] = edges[j] ;
         edges[j] = wid      ;
         j-- ;
      } else {
         w->mark = 'X' ;
         i++ ;
      }
   }
   v->nadj = j + 1 ;
} else {
/*
   ----------------------------------------------------
   v is not a leaf, merge its subtrees' boundaries
   with its uncovered edge list to get the new boundary
   ----------------------------------------------------
*/
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n vertex %d is not a leaf", v->id) ;
      fprintf(info->msgFile, "\n  vertex %d, subtrees :", v->id) ;
      IP_fp80(info->msgFile, v->subtrees, 20) ;
      fflush(info->msgFile) ;
   }
   v->status = 'E' ;
   nadj = 0 ;
   while ( (ip = v->subtrees) != NULL ) {
      if ( info->msglvl > 3 ) {
         fprintf(info->msgFile, "\n    subtree %d, ip(%p)<%d,%p>",
                 ip->val, ip, ip->val, ip->next) ;
         fflush(info->msgFile) ;
      }
      uid    = ip->val ;
      u      = msmd->vertices + uid ;
      u->par = v ;
      nbnd   = u->nadj ;
      bnd    = u->adj  ;
      if ( info->msglvl > 3 ) {
         fprintf(info->msgFile, "\n    bnd of adj subtree %d :", u->id) ;
         IVfp80(info->msgFile, nbnd, bnd, 25, &ierr) ;
         fflush(info->msgFile) ;
      }
      for ( i = 0 ; i < nbnd ; i++ ) {
         wid = bnd[i] ;
         w   = msmd->vertices + wid ;
         if ( w->mark == 'O' && w->status != 'I' ) {
            w->mark = 'X' ;
            adj[nadj++] = wid ;
         }
      }
      if ( u->status == 'E' ) {
/*
         ------------------------------------------
         u is not a leaf, free its boundary storage
         ------------------------------------------
*/
         IVfree(u->adj) ;
         info->nbytes -= u->nadj * sizeof(int) ;
      }
      u->adj  = NULL ;
      u->nadj =   0  ;
/*
      --------------------------------------
      put this IP structure on the free list
      --------------------------------------
*/
      v->subtrees = ip->next    ;
      ip->val     = -1          ;
      ip->next    = msmd->freeIP ;
      msmd->freeIP = ip          ;
      if ( info->msglvl > 3 ) {
         fprintf(info->msgFile, 
                 "\n   v->subtrees = %p, msmd->freeIP = %p",
                 v->subtrees, msmd->freeIP) ;
         fflush(info->msgFile) ;
      }
   }
/*
   ------------------------------------------------
   merge all uncovered edges into the boundary list
   ------------------------------------------------
*/
   nedge = v->nadj ;
   edges = v->adj  ;
   for ( i = 0 ; i < nedge ; i++ ) {
      wid = edges[i] ;
      w   = msmd->vertices + wid ;
      if ( w->mark == 'O' && w->status != 'I' ) {
         w->mark     = 'X' ;
         adj[nadj++] = wid ;
      }
   }
/*
   ----------------------------------------------------------
   if boundary is not empty, allocate new storage for entries
   ----------------------------------------------------------
*/
   v->nadj = nadj ;
   if ( nadj > 0 ) {
      v->adj = IVinit(nadj, -1) ;
      IVcopy(nadj, v->adj, adj) ;
      info->nbytes += nadj*sizeof(int) ;
      if ( info->maxnbytes < info->nbytes ) {
         info->maxnbytes = info->nbytes ;
      }
   } else {
      v->adj = NULL ;
   }
}
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n    bnd(%d) :", v->id) ;
   if ( v->nadj > 0 ) {
      IVfp80(info->msgFile, v->nadj, v->adj, 17, &ierr) ;
   }
   fflush(info->msgFile) ;
}
/*
   ----------------------------------------------
   for each boundary vertex
      1. add v to subtree list
      2. put v on reach set if not already there
      3. unmark and add weight to boundary weight
   ----------------------------------------------
*/
nbnd = v->nadj ;
bnd  = v->adj  ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n %d's bnd :", v->id) ;
   IVfp80(info->msgFile, nbnd, bnd, 12, &ierr) ;
   fflush(info->msgFile) ;
}
wght = 0 ;
for ( i = 0 ; i < nbnd ; i++ ) {
   wid = bnd[i] ;
   w   = msmd->vertices + wid ;
   if ( info->msglvl > 4 ) {
      fprintf(info->msgFile, "\n   adjacent vertex %d", w->id) ;
      fflush(info->msgFile) ;
   }
/*
   -------------------------------
   add v to the subtree list for w
   -------------------------------
*/
   if ( (ip = msmd->freeIP) == NULL ) {
      if ( info->msglvl > 2 ) {
         fprintf(info->msgFile, "\n   need to get more IP objects") ;
         fflush(info->msgFile) ;
      }
/*
      -------------------------------------------------
      no more free IP structures, allocate more storage
      -------------------------------------------------
*/
      if ( (ip = IP_init(msmd->incrIP, IP_FORWARD)) == NULL ) {
         fprintf(stderr, "\n fatal error in MSMD_eliminateVtx%p,%p,%p)"
                 "\n unable to allocate more IP objects",
                 msmd, v, info) ;
         exit(-1) ;
      }
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, "\n   old baseIP = %p", msmd->baseIP) ;
         fprintf(info->msgFile, "\n   new baseIP = %p", ip) ;
         fflush(info->msgFile) ;
      }
      ip->next = msmd->baseIP ;
      msmd->baseIP = ip ;
      info->nbytes += msmd->incrIP*sizeof(struct _IP) ;
      if ( info->maxnbytes < info->nbytes ) {
         info->maxnbytes = info->nbytes ;
      }
      ip = msmd->freeIP = msmd->baseIP + 1 ;
      if ( info->msglvl > 2 ) {
         fprintf(info->msgFile, "\n   all set") ;
         fflush(info->msgFile) ;
      }
   }
   msmd->freeIP = ip->next ;
   ip->val     = v->id    ;
   ip->next    = NULL     ;
   for ( ip2 = w->subtrees, prev = NULL ; 
         ip2 != NULL && ip2->val > ip->val ;
         ip2 = ip2->next ) {
      prev = ip2 ;
   }
   if ( prev == NULL ) {
      w->subtrees = ip ;
   } else {
      prev->next = ip ;
   }
   ip->next = ip2 ;
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n %d's subtrees :", w->id) ;
      IP_fp80(info->msgFile, w->subtrees, 15) ;
      fflush(info->msgFile) ;
   }
/*
   --------------------------------
   add w to reach list if necessary
   --------------------------------
*/
   if ( info->msglvl > 4 ) {
      fprintf(info->msgFile, "\n    status[%d] = %c", wid, w->status) ;
      fflush(info->msgFile) ;
   }
   switch ( w->status ) {
   case 'D' :
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, ", remove from heap") ;
         fflush(info->msgFile) ;
      }
      IIheap_remove(msmd->heap, wid) ;
   case 'O' :
   case 'B' :
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, ", add to reach set, nreach = %d",
                 IV_size(reachIV) + 1) ;
         fflush(info->msgFile) ;
      }
      IV_push(reachIV, wid) ;
      w->status = 'R' ;
   case 'R' :
      break ;
   case 'I' :
      break ;
   default :
      fprintf(stderr, "\n error in MSMD_eliminateVtx(%p,%p,%p)"
              "\n status[%d] = '%c'\n",
              msmd, v, info, wid, w->status) ;
      fprintf(stderr, "\n msmd->nvtx = %d", msmd->nvtx) ;
      exit(-1) ;
   }
/*
   --------------------------------
   unmark the boundary vertices and
   store the weight of the boundary
   --------------------------------
*/
   w->mark = 'O' ;
   wght += w->wght ;
}
/*
   ------------------------------------
   unmark v and set its boundary weight
   ------------------------------------
*/
v->mark    = 'O'  ;
v->bndwght = wght ;

return ; }

/*--------------------------------------------------------------------*/
