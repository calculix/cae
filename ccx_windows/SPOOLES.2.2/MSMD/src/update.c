/*  update.c  */

#include "../MSMD.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   purpose -- to update vertices in the reach set

   created -- 96feb25, cca
   ----------------------------------------------
*/
void
MSMD_update ( 
   MSMD       *msmd,
   MSMDinfo   *info
) {
int       ii, jj, nreach, vid, wght ;
int       *reach ;
IP        *ip ;
IV        *reachIV ;
MSMDvtx   *v ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_update(%p,%p)"
           "\n bad input\n", msmd, info) ;
   exit(-1) ;
}
if ( info->msglvl > 4 ) {
   fprintf(info->msgFile, 
           "\n inside MSMD_update(%p,%p), nreach = %d", 
           msmd, info, IV_size(&msmd->reachIV)) ;
   fflush(info->msgFile) ;
}
/*
   ---------------------------------
   if the reach set is empty, return
   ---------------------------------
*/
reachIV = &msmd->reachIV ;
if ( (nreach = IV_size(reachIV)) == 0 ) {
   return ; 
}
reach = IV_entries(reachIV) ;
if ( info->msglvl > 4 ) {
   for ( ii = 0 ; ii < nreach ; ii++ ) {
      vid = reach[ii] ;
      v   = msmd->vertices + vid ;
      MSMDvtx_print(v, info->msgFile) ;
   }
   fflush(info->msgFile) ;
}
/*
   ----------------------------
   switch over the update types
   ----------------------------
*/
if ( info->prioType == 2 ) {
/*
   -------------------
   approximate updates
   -------------------
*/
   for ( ii = 0 ; ii < nreach ; ii++ ) {
      vid = reach[ii] ;
      v   = msmd->vertices + vid ;
      if ( v->status == 'I' ) {
/*
         -------------------------
         node is indistinguishable
         -------------------------
*/
         continue ;
      } else if ( v->status == 'R' ) {
/*
         ---------------------------------------------
         internal node needs approximate degree update
         ---------------------------------------------
*/
         wght = MSMD_approxDegree(msmd, v, info) ;
         if ( info->msglvl > 3 ) {
            fprintf(info->msgFile, 
                    "\n inserting %d with priority %d into heap",
                    vid, wght) ;
            fflush(info->msgFile) ;
         }
         IIheap_insert(msmd->heap, vid, wght) ;
         v->status = 'D' ;
      }
   }
} else if ( info->prioType == 0 ) {
/*
   --------------------------------------------------------
   maximal independent set elimination, set degrees to zero
   --------------------------------------------------------
*/
   for ( ii = 0 ; ii < nreach ; ii++ ) {
      vid = reach[ii] ;
      v   = msmd->vertices + vid ;
      if ( v->status != 'I' ) {
         IIheap_insert(msmd->heap, vid, 0) ;
         v->status = 'D' ;
      }
   }
} else {
/*
   ---------------------------------------------------------
   exact or half and half updates, process 2-adj nodes first
   ---------------------------------------------------------
*/
   for ( ii = 0, jj = 0 ; ii < nreach ; ii++ ) {
      vid = reach[ii] ;
      v   = msmd->vertices + vid ;
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, "\n v = %d, stage = %d, status = %c",
                 v->id, v->stage, v->status) ;
         fflush(info->msgFile) ;
      }
      if ( v->status == 'R' ) {
/*
         ---------------------------------
         internal node needs degree update
         ---------------------------------
*/
         if (  v->nadj == 0
            && (ip = v->subtrees) != NULL
            && (ip = ip->next) != NULL
            && ip->next == NULL ) {
/*
            --------------------------------------------------------
            v is adjacent to two subtrees and has no uncovered edges
            --------------------------------------------------------
*/
            if ( info->msglvl > 4 ) {
               fprintf(info->msgFile, ", 2-adj vertex") ;
               fflush(info->msgFile) ;
            }
            wght = MSMD_exactDegree2(msmd, v, info) ;
            if ( info->msglvl > 4 ) {
               fprintf(info->msgFile, 
                  "\n   2-adj, inserting %d with priority %d into heap",
                  vid, wght) ;
               fflush(info->msgFile) ;
            }
            IIheap_insert(msmd->heap, vid, wght) ;
            v->status = 'D' ;
         } else {
/*
            -----------------------------------------
            v is adjacent two or more subtrees and/or 
            has uncovered edges, keep in reach list
            -----------------------------------------
*/
            reach[jj++] = vid ;
         }
      }
   }
   nreach = jj ;
/*
   ----------------------------------------------------
   second pass, update the remaining nodes in the reach
   ----------------------------------------------------
*/
   for ( ii = 0 ; ii < nreach ; ii++ ) {
      vid = reach[ii] ;
      v   = msmd->vertices + vid ;
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, "\n v = %d, stage = %d, status = %c",
                 v->id, v->stage, v->status) ;
         fflush(info->msgFile) ;
      }
      if ( v->status == 'R' ) {
         if ( info->prioType == 1 ) {
            wght = MSMD_exactDegree3(msmd, v, info) ;
         } else {
            wght = MSMD_approxDegree(msmd, v, info) ;
         }
         if ( info->msglvl > 4 ) {
            fprintf(info->msgFile, 
                 "\n   3-adj, inserting %d with priority %d into heap",
                    vid, wght) ;
            fflush(info->msgFile) ;
         }
         IIheap_insert(msmd->heap, vid, wght) ;
         v->status = 'D' ;
      }
   }
}
/*
msmd->nreach = 0 ;
*/
IV_setSize(reachIV, nreach) ;
/*
   ------------------------------------
   optionally print out the degree heap
   ------------------------------------
*/
if ( info->msglvl > 4 ) {
   fprintf(info->msgFile, "\n degree heap") ;
   IIheap_print(msmd->heap, info->msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to compute the exact boundary size of a node 
              adjacent to only two eliminated vertices

   created -- 96feb25, cca
   -------------------------------------------------------
*/
int
MSMD_exactDegree2 ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) {
int       bndwght, i, j, usize0, usize1, wid ;
int       *uadj0, *uadj1 ;
MSMDvtx   *u0, *u1, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_exactDegree2(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
if ( v->subtrees == NULL ) {
   fprintf(stderr, "\n\n 1. error in MSMD_exactDegree2(%p,%p,%p)"
           "\n v->subtrees == NULL\n", msmd, v, info) ;
   exit(-1) ;
}
if ( v->subtrees->next == NULL ) {
   fprintf(stderr, "\n\n 1. error in MSMD_exactDegree2(%p,%p,%p)"
           "\n v->subtrees->next == NULL\n", msmd, v, info) ;
   exit(-1) ;
}
/*
   -------------------------------------------------------
   this vertex is adjacent to only two eliminated vertices
   -------------------------------------------------------
*/
u0     = msmd->vertices + v->subtrees->val ;
usize0 = u0->nadj ; 
uadj0  = u0->adj  ;
if ( usize0 == 0 || uadj0 == NULL ) {
   fprintf(stderr, "\n\n 1. error in MSMD_exactDegree2(%p,%p,%p)"
           "\n usize0 = %d, uadj0 = %p"
           "\n bad adjacency list for %d\n ", 
           msmd, v, info, usize0, uadj0, u0->id) ;
   MSMDvtx_print(u0, stderr) ;
   exit(-1) ;
}
u1     = msmd->vertices + v->subtrees->next->val ;
usize1 = u1->nadj ; 
uadj1  = u1->adj  ;
if ( usize1 == 0 || uadj1 == NULL ) {
   fprintf(stderr, "\n\n 2. error in MSMD_exactDegree2(%p,%p,%p)"
           "\n usize1 = %d, uadj1 = %p"
           "\n bad adjacency list for %d\n ", 
           msmd, v, info, usize1, uadj1, u1->id) ;
   MSMDvtx_print(u1, stderr) ;
   exit(-1) ;
}
/*
   -----------------------------------------------
   mark all nodes adjacent to the first eliminated 
   vertex and add their size to the boundary. 
   -----------------------------------------------
*/
v->mark = 'X' ;
bndwght = 0 ;
i = 0, j = usize0 - 1 ;
while ( i <= j ) {
   wid = uadj0[i] ;
   w   = msmd->vertices + wid ;
   if ( w->status == 'I' ) {
      uadj0[i] = uadj0[j] ;
      uadj0[j] = wid ;
      j-- ;
   } else {
      if ( w->mark != 'X' ) {
         w->mark = 'X' ;
         bndwght += w->wght ;
         if ( info->msglvl > 5 ) {
            fprintf(info->msgFile, 
                    "\n    %d : adding %d with weight %d to boundary",
                    u0->id, w->id, w->wght) ;
            fflush(info->msgFile) ;
         }
      }
      i++ ;
   }
}
u0->nadj = j + 1 ;
/*
   ---------------------------------------------------
   loop over all nodes adjacent to the second 
   eliminated vertex. if a node is marked, it is 
   outmatched by v, else add its size to the boundary.
   ---------------------------------------------------
*/
i = 0, j = usize1 - 1 ;
while ( i <= j ) {
   wid = uadj1[i] ;
   w   = msmd->vertices + wid ;
   if ( w == v ) {
      i++ ;
   } else {
      if ( w->status == 'I' ) {
         uadj1[i] = uadj1[j] ;
         uadj1[j] = wid ;
         j-- ;
      } else {
         if ( w->mark != 'X' ) {
            bndwght += w->wght ;
            if ( info->msglvl > 5 ) {
               fprintf(info->msgFile, 
                      "\n    %d : adding %d with weight %d to boundary",
                       u1->id, w->id, w->wght) ;
               fflush(info->msgFile) ;
            }
         } else if ( w->status == 'R' ) {
            if (  info->msglvl > 2 ) {
               fprintf(info->msgFile, 
                       "\n    %6d is outmatched by %6d", w->id, v->id) ;
               fflush(info->msgFile) ;
            }
            w->status = 'O' ;
            info->stageInfo->noutmtch++ ;
         }
         i++ ;
      }
   }
}
u1->nadj = j + 1 ;
/*
   -------------------------------------------
   unmark the nodes in the first adjacency set
   -------------------------------------------
*/
usize0  = u0->nadj ;
v->mark = 'O' ;
for ( j = 0 ; j < usize0 ; j++ ) {
   wid = uadj0[j] ;
   w   = msmd->vertices + wid ;
   w->mark = 'O' ;
}
/*
   -------------------------------------
   increment the number of exact updates
   -------------------------------------
*/
info->stageInfo->nexact2++ ;

return(bndwght) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to compute the exact boundary size of a node that
              is not adjacent to only two eliminated vertices

   created -- 96feb25, cca
   ------------------------------------------------------------
*/
int
MSMD_exactDegree3 ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) {
int       bndwght, i, ierr, j, nadj, uid, usize, vsize, wid ;
int       *adj, *uadj, *vadj ;
IP        *ip ;
MSMDvtx   *u, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_exactDegree3(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
adj = IV_entries(&msmd->ivtmpIV) ;
/*
   ---------------------------------
   accumulate the boundary set for v
   ---------------------------------
*/
v->mark = 'X' ;
nadj = 0 ;
for ( ip = v->subtrees ; ip != NULL ; ip = ip->next ) {
   uid = ip->val ;
   u = msmd->vertices + uid ;
   usize = u->nadj ;
   uadj  = u->adj  ;
   i = 0, j = usize - 1 ;
   while ( i <= j ) {
      wid = uadj[i] ;
      w   = msmd->vertices + wid ;
      if ( w->status == 'I' ) {
         uadj[i] = uadj[j] ;
         uadj[j] = wid     ;
         j-- ;
      } else {
         i++ ;
         if ( w->mark != 'X' ) {
            w->mark = 'X' ;
            adj[nadj++] = wid ;
         }
      }
   }
   u->nadj = j + 1 ;
}
vsize = v->nadj ;
vadj  = v->adj  ;
i = 0, j = vsize - 1 ;
while ( i <= j ) {
   uid = vadj[i] ;
   u = msmd->vertices + uid ;
   if ( u->status == 'I' ) {
      vadj[i] = vadj[j] ;
      vadj[j] = uid ;
      j-- ;
   } else {
      if ( u->mark != 'X' ) {
         u->mark = 'X' ;
         adj[nadj++] = uid ;
      }
      i++ ;
   }
}
v->nadj = j + 1 ;
if ( info->msglvl > 4 ) {
   fprintf(info->msgFile, "\n vadj(%d) :", v->id) ;
   IVfp80(info->msgFile, vsize, vadj, 12, &ierr) ;
   fflush(info->msgFile) ;
}
/*
   ---------------------------------------------------
   the boundary of v is found in the list adj, compute 
   the boundary size and unmark the boundary vertices
   ---------------------------------------------------
*/
bndwght = 0 ;
for ( i = 0 ; i < nadj ; i++ ) {
   uid = adj[i] ;
   u   = msmd->vertices + uid ;
   bndwght += u->wght ;
   u->mark = 'O' ;
}
v->mark = 'O' ;
/*
   -------------------------------------
   increment the number of exact updates
   -------------------------------------
*/
info->stageInfo->nexact3++ ;

return(bndwght) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to compute the approximate degree of a vertex

   created -- 96feb25, cca
   --------------------------------------------------------
*/
int
MSMD_approxDegree ( 
   MSMD       *msmd,
   MSMDvtx    *v,
   MSMDinfo   *info
) {
int       bndwght, i, uid, vsize ;
int       *vadj ;
IP        *ip ;
MSMDvtx   *u ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || v == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_approxDegree(%p,%p,%p)"
           "\n bad input\n", msmd, v, info) ;
   exit(-1) ;
}
/*
   ------------------------------------------------
   accumulate the approximate boundary weight for v
   ------------------------------------------------
*/
bndwght = 0 ;
for ( ip = v->subtrees ; ip != NULL ; ip = ip->next ) {
   bndwght += msmd->vertices[ip->val].bndwght - v->wght ;
}
vsize = v->nadj ;
vadj  = v->adj  ;
for ( i = 0 ; i < vsize ; i++ ) {
   uid = vadj[i] ;
   u = msmd->vertices + uid ;
   if ( u != v && u->status != 'I' ) {
      bndwght += u->wght ;
   }
}
/*
   -------------------------------------------
   increment the number of approximate updates
   -------------------------------------------
*/
info->stageInfo->napprox++ ;

return(bndwght) ; }

/*--------------------------------------------------------------------*/
