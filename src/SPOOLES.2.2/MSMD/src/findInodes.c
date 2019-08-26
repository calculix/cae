/*  findInodes2.c  */

#include "../MSMD.h"

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- to find indistinguishable nodes in the reach set

   flag = 0 --> return
   flag = 1 --> check out nodes that are 2-adj
   flag = 2 --> check out nodes that are both 2-adj and not

   note: the reach set is not changed.

   created -- 96feb15, cca
   modified -- 97feb07, cca
      very tricky "bug"
      was : sum += ip->val for subtrees
            sum += IVsum(nvedge, vedges) for uncovered edges
      now : sum += ip->val + 1 for subtrees
            sum += IVsum(nvedge, vedges) + nvedge for uncovered edges
      checksums were "wrong" due to vertex 0 adding nothing
      to the checksum. beware 0-indexing.
   ------------------------------------------------------------------
*/
void
MSMD_findInodes ( 
   MSMD       *msmd,
   MSMDinfo   *info
) {
int             first, flag, i, ierr, iv, iw, j, k, keepon, nlist, 
                nreach, nvedge, sum, vid, vchk, vcount, wid ;
int             *chk, *list, *reach, *vedges, *wedges ;
IP              *ip, *ipv, *ipw, *vsubtrees ;
MSMDvtx         *v, *w ;
/*
   ---------------
   check the input
   ---------------
*/
if ( msmd == NULL || info == NULL ) {
   fprintf(stderr, "\n fatal error in MSMD_findInodes(%p,%p)"
           "\n bad input\n", msmd, info) ;
   exit(-1) ;
}
if ( (flag = info->compressFlag % 4) == 0 ) {
/*
   ---------------------------------------
   no compression requested, simple return
   ---------------------------------------
*/
   return ;
}
/*
   ---------------------------------
   if the reach set is empty, return
   ---------------------------------
*/
if ( (nreach = IV_size(&msmd->reachIV)) == 0 ) {
   return ; 
}
/*
reach = msmd->reach ;
*/
reach = IV_entries(&msmd->reachIV) ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n inside MSMD_findInodes(%p)"
           "\n reach(%d) :", msmd, nreach) ;
   IVfp80(info->msgFile, nreach, reach, 10, &ierr);
   fflush(info->msgFile) ;
}
/*
   -------------------------------------------------------
   load the front of the reach set with nodes to be tested
   -------------------------------------------------------
*/
chk = IV_entries(&msmd->ivtmpIV) ;
list = reach ;
if ( flag == 1 ) {
/*
   -------------------------------------------
   work only with nodes adjacent to 2 subtrees
   -------------------------------------------
*/
   i = 0 ; j = nreach - 1 ;
   while ( i <= j ) {
      vid = list[i] ;
      v = msmd->vertices + vid ;
      if (  v->nadj != 0 || (ip = v->subtrees) == NULL
        || (ip = ip->next) == NULL || ip->next != NULL ) {
/*
         --------------------------------
         vertex is not 2-adj, swap to end
         --------------------------------
*/
         list[i] = list[j] ;
         list[j] = vid     ;
         j-- ;
      } else {
/*
         --------------------------
         vertex is 2-adj, keep here
         --------------------------
*/
         i++ ;
      }
   }
   nlist = j + 1 ;
} else {
/*
   ---------------------------------
   put all reached nodes in the list
   ---------------------------------
*/
   nlist = nreach ;
}
if ( nlist == 0 ) {
   return ; 
}
/*
   -----------------------------------------------------
   compute the the checksums and count adjacent subtrees
   for all vertices in the list
   -----------------------------------------------------
*/
for ( k = 0 ; k < nlist ; k++ ) {
   vid    = list[k] ;
   v      = msmd->vertices + vid ;
   vcount = 0 ;
   sum    = 0 ;
   if ( info->msglvl > 4 ) {
      fprintf(info->msgFile, "\n vertex %d", vid) ;
      fflush(info->msgFile) ;
   }
   for ( ipv = v->subtrees ; ipv != NULL ; ipv = ipv->next ) {
/*
      ------------------------------------
      add adjacent subtree to the checksum
      ------------------------------------
*/
      sum += ipv->val + 1 ;
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, "\n    adjacent subtree %d, sum = %d",
                 ipv->val, sum) ;
         fflush(info->msgFile) ;
      }
      vcount++ ;
   }
   if ( (nvedge = v->nadj) > 0 && (vedges = v->adj) != NULL ) {
      sum += IVsum(nvedge, vedges) + nvedge ;
      if ( info->msglvl > 4 ) {
         fprintf(info->msgFile, "\n    %d adjacent edges :",
                 nvedge) ;
         IVfp80(info->msgFile, nvedge, vedges, 20, &ierr) ;
         fprintf(info->msgFile, " : sum = %d", sum) ;
         fflush(info->msgFile) ;
      }
      IVqsortUp(nvedge, vedges) ;
   }
/*
   -----------------
   save the checksum
   -----------------
*/
   chk[k] = sum ;
}
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n before sort, list array") ;
   fflush(info->msgFile) ;
   IVfp80(info->msgFile, nlist, list, 80, &ierr) ;
   fflush(info->msgFile) ;
   fprintf(info->msgFile, "\n chk array") ;
   fflush(info->msgFile) ;
   IVfp80(info->msgFile, nlist, chk, 80, &ierr) ;
   fflush(info->msgFile) ;
}
/*
   -----------------------------------------------------
   sort the vertices in the reach set by their checksums
   -----------------------------------------------------
*/
IV2qsortUp(nlist, chk, list) ;
if ( info->msglvl > 3 ) {
   fprintf(info->msgFile, "\n after sort, reach array") ;
   IVfp80(info->msgFile, nlist, list, 80, &ierr) ;
   fprintf(info->msgFile, "\n chk array") ;
   IVfp80(info->msgFile, nlist, chk, 80, &ierr) ;
   fflush(info->msgFile) ;
}
/*
   ----------------------------------------
   detect and purge indistinguishable nodes
   ----------------------------------------
*/
for ( iv = 0 ; iv < nlist ; iv++ ) {
   vid = list[iv] ;
   v   = msmd->vertices + vid ;
   if ( v->status == 'I' ) {
/*
      -----------------------------------------------------
      vertex has been found indistinguishable, skip to next
      -----------------------------------------------------
*/
      continue ;
   }
/*
   ---------------------------
   test against other vertices
   ---------------------------
*/
   vchk      = chk[iv]     ;
   nvedge    = v->nadj     ;
   vedges    = v->adj      ;
   vsubtrees = v->subtrees ;
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, 
              "\n checking out v = %d, vchk = %d, status = %c",
              v->id, vchk, v->status) ;
      fflush(info->msgFile) ;
   }
/*
   ---------------------------------------------------
   check v against all vertices with the same checksum
   ---------------------------------------------------
*/
   if ( info->msglvl > 3 ) {
      fprintf(info->msgFile, "\n checking out v = %d, status = %d",
              v->id, v->stage) ;
      fflush(info->msgFile) ;
   }
   first =  1 ;
   for ( iw = iv + 1 ; iw < nlist && chk[iw] == vchk ; iw++ ) {
      wid = reach[iw] ;
      w   = msmd->vertices + wid ;
      if ( info->msglvl > 3 ) {
         fprintf(info->msgFile, 
                 "\n     w = %d, status = %c, status = %d",
                 w->id, w->status, w->stage) ;
         fflush(info->msgFile) ;
      }
      if (  w->status == 'I' 
         || v->stage != w->stage
         || nvedge != w->nadj ) {
/*
         ------------------------------------
         w has been found indistinguishable
            or
         v and w do not lie on the same stage
            or
         edge counts are not the same
         ------------------------------------
*/
         continue ;
      }
/*
      ----------------------------------------
      w and v check out so far, check to see 
      if all vertices adjacent to w are marked
      ----------------------------------------
*/
      if ( info->msglvl > 3 ) {
         fprintf(info->msgFile, 
                 "\n    checking %d against %d", wid, vid) ;
         fflush(info->msgFile) ;
      }
/*
      ---------------------------------------------------------------
      check to see if the subtree lists and edge lists are indentical
      ---------------------------------------------------------------
*/
      info->stageInfo->ncheck++ ;
      keepon = 1 ;
      ipv = vsubtrees ;
      ipw = w->subtrees ;
      while ( ipv != NULL && ipw != NULL ) {
         if ( ipv->val != ipw->val ) {
            keepon = 0 ;
            break ;
         }
         ipv = ipv->next ;
         ipw = ipw->next ;
      }
      if ( keepon == 1 ) {
         wedges = w->adj ;
         for ( k = 0 ; k < nvedge ; k++ ) {
            if ( vedges[k] != wedges[k] ) {
               keepon = 0 ;
               break ;
            }
         }
      }
      if ( keepon == 1 ) {
/*
         ---------------------------------------------
         w and v are indistinguishable, merge w into v
         ---------------------------------------------
*/
         if ( info->msglvl > 1 ) {
            fprintf(info->msgFile, 
                    "\n %d absorbs %d, wght = %d, status[%d] = %c", 
                    v->id, w->id, w->wght, w->id, w->status) ;
            fflush(info->msgFile) ;
         }
         v->wght   += w->wght  ;
         w->wght   =  0 ;
         w->status =  'I' ;
         w->nadj   =  0 ;
         w->adj    = NULL ;
         w->par    =  v ;
         if ( (ipw = w->subtrees) != NULL ) {
            while ( ipw->next != NULL ) {
               ipw = ipw->next ;
            }
            ipw->next   = msmd->freeIP ;
            msmd->freeIP = ipw ;
            w->subtrees = NULL ;
         }
         info->stageInfo->nindst++ ;
      }
   }
}
if ( info->msglvl > 4 ) {
   fprintf(info->msgFile, 
         "\n MSMD_findInodes(%p), all done checking the nodes", msmd) ;
   fflush(info->msgFile) ;
}

return ; }

/*--------------------------------------------------------------------*/
