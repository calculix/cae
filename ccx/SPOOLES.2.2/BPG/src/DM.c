/*  DM.c  */

#include "../../BPG.h"


/*--------------------------------------------------------------------*/
/*
   -----------------
   static prototypes
   -----------------
*/
static void
nonunitDM ( BPG *bpg, int dmflags[], int stats[], 
            int msglvl, FILE *msgFile) ;
static void
unitDM ( BPG *bpg, int dmflags[], int stats[], 
         int msglvl, FILE *msgFile) ;
static int
nonunitFindExposedNode ( BPG *bpg, int uexp, int list[], int link[],
                         int mark[], int tag, int nvexp[], IVL *ivl,
                         int msglvl, FILE *msgFile ) ;
static int
nonunitFindNmatch ( BPG *bpg, int uexp, int vexp, int nvexp[],
                    int link[], IVL *ivl, int msglvl, FILE *msgFile ) ;
static void
nonunitFlipEdges ( BPG *bpg, int uexp, int vexp, int nmatch,
                   int nvexp[], int link[], IVL *ivl, int msglvl,
                   FILE *msgFile ) ;
static void
nonunitSetFlags ( BPG *bpg, int list[], int nvexp[], int mark[],
                  IVL *ivl, int dmflags[], int stats[], int msglvl,
                  FILE *msgFile ) ;
static void
unitFindMaxMatch ( BPG *bpg, int mate[], int msglvl, FILE *msgFile ) ;
static int
unitAugmentingPath ( BPG *bpg, int uexp, int mate[], int mark[],
                   int link[], int list[], int msglvl, FILE *msgFile ) ;
static void
unitSetFlags ( BPG *bpg, int mate[], int dmflags[], int stats[],
               int msglvl, FILE *msgFile ) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   compute the generalized Dulmadge-Mendolsohn decomposition

   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 95oct12, cca
   ---------------------------------------------------------
*/
void
BPG_DMdecomposition (
   BPG   *bpg,
   int   dmflags[],
   int   stats[],
   int   msglvl,
   FILE   *msgFile
) {
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n inside BPG_DMdecomposition") ;
   fflush(msgFile) ;
}
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || dmflags == NULL || stats == NULL 
   || (msglvl > 0 && msgFile == NULL) ) {
   fprintf(stderr, 
           "\n fatal error in BPG_DMdecomposition(%p,%p,%p,%d,%p)"
           "\n bad input\n", bpg, dmflags, stats, msglvl, msgFile) ;
   exit(-1) ;
}
if ( bpg->graph->type % 2 == 0 ) {
/*
   -----------------
   unit weight graph 
   -----------------
*/
   unitDM(bpg, dmflags, stats, msglvl, msgFile) ;
} else {
/*
   --------------------
   nonunit weight graph 
   --------------------
*/
   nonunitDM(bpg, dmflags, stats, msglvl, msgFile) ;
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n leaving BPG_DMdecomposition") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   compute the generalized Dulmadge-Mendolsohn decomposition
   for a nonunit weight graph.

   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 95oct12, cca
   ---------------------------------------------------------
*/
static void
nonunitDM (
   BPG   *bpg,
   int   dmflags[],
   int   stats[],
   int   msglvl,
   FILE   *msgFile
) {
int   maxnXnY, nmatch, nX, nY, tag, x, xexp, y, yexp ;
int   *link, *list, *mark, *nvexp ;
IVL   *ivl ;

if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n inside nonunitDM") ;
   fflush(msgFile) ;
}
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || dmflags == NULL || stats == NULL ) {
   fprintf(stderr, "\n fatal error in nonunitDM(%p,%p,%p)"
           "\n bad input\n", bpg, dmflags, stats) ;
   exit(-1) ;
}
nX = bpg->nX ;
nY = bpg->nY ;
maxnXnY = ( nX >= nY ) ? nX : nY ;
/*
   -----------------------
   create the working data
   -----------------------
*/
list  = IVinit(maxnXnY, -1) ;
link  = IVinit(nX+nY,   -1) ;
mark  = IVinit(nX+nY,   -1) ;
nvexp = IVinit(nX+nY,    0) ;
IVcopy(nX + nY, nvexp, bpg->graph->vwghts) ;
ivl = IVL_new() ;
IVL_setDefaultFields(ivl) ;
IVL_init3(ivl, IVL_CHUNKED, nX + nY, nvexp) ;
/*
   --------
   DEBUG!!!
   --------
*/
/*
{ int  v ;
for ( v = 0 ; v < nX + nY ; v++ ) {
   IVshuffle(bpg->g->adjIVL->sizes[v], bpg->g->adjIVL->p_vec[v], 57) ;
}
}
*/
/*
   ------------------------
   loop over the x vertices
   ------------------------
*/
for ( xexp = 0, tag = 1 ; xexp < nX ; xexp++ ) {
   while ( nvexp[xexp] > 0 ) {
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n\n checking out %d, nvexp[%d] = %d",
                 xexp, xexp, nvexp[xexp]) ;
      }
      yexp = nonunitFindExposedNode(bpg, xexp, list, link, mark,
                                    tag, nvexp, ivl, msglvl, msgFile) ;
      tag++ ;
      if ( yexp == -1 ) {
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, 
                 "\n    exposed node %d, no alternating path", xexp) ;
            fflush(msgFile) ;
         }
         break ;
      } else {
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, 
                 "\n    exposed node %d, alternating node ends at %d",
                 xexp, yexp) ;
            fprintf(msgFile, "\n    path : %d", yexp) ;
            x = link[yexp] ; 
            while ( x != xexp ) {
               y = link[x] ;
               fprintf(msgFile, " --> %d ==> %d", x, y) ;
               x = link[y] ;
            }
            fprintf(msgFile, " --> %d", xexp) ;
            fflush(msgFile) ;
         }
         nmatch = nonunitFindNmatch(bpg, xexp, yexp, nvexp, link, 
                                    ivl, msglvl, msgFile) ; 
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n nmatch = %d", nmatch) ;
         }
         nonunitFlipEdges(bpg, xexp, yexp, nmatch, nvexp, link, 
                          ivl, msglvl, msgFile) ; 
      }
   }
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n match IVL") ;
   IVL_writeForHumanEye(ivl, msgFile) ;
}
/*
   -------------
   set the flags
   -------------
*/
nonunitSetFlags(bpg, list, nvexp, mark, ivl, 
                dmflags, stats, msglvl, msgFile) ; 
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list)  ;
IVfree(link)  ;
IVfree(mark)  ;
IVfree(nvexp) ;
IVL_free(ivl) ;

if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n leaving nonunitDM") ;
   fflush(msgFile) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------------
   find exposed node that is the endpoint of 
   an alternating path that starts with uexp
   this method uses a breadth-first traversal

   uexp    -- exposed node to be the start of the alternating path
   list    -- list vector, size max(nX, nY)
   link    -- link vector, used to define alternating path, size nX + nY
   mark    -- mark vector, used to keep track of nodes in the tree
   tag     -- unique tag for this call
   nvexp   -- nvexp[u] = # of exposed vertices in u
   ivl     -- IVL object to hold matched edges

   return value -- exposed node if exists, otherwise -1

   created -- 95oct14, cca
   ---------------------------------------------------------------------
*/
static int
nonunitFindExposedNode (
   BPG   *bpg,
   int   uexp,
   int   list[],
   int   link[],
   int   mark[],
   int   tag,
   int   nvexp[],
   IVL   *ivl,
   int   msglvl,
   FILE   *msgFile
) {
int   ierr, ii, jj, last, now, u, unew, usize, v, vsize ;
int   *uadj, *vind ;

if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n inside nonunitFindExposedNode(%d)", uexp) ;
   fflush(msgFile) ;
}
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || uexp < 0 || uexp > (bpg->nX + bpg->nY)
   || list == NULL || link == NULL 
   || mark == NULL || nvexp == NULL || ivl == NULL ) {
   fprintf(stderr, 
    "\n fatal error in nonunitFindExposedNode2(%p,%d,%p,%p,%p,%d,%p,%p)"
    "\n bad input\n",
    bpg, uexp, list, link, mark, tag, nvexp, ivl) ;
   exit(-1) ;
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n\n working on exposed %d, nvexp %d", 
           uexp, nvexp[uexp]) ;
}
/*
   -----------------------------------
   load the list with the exposed node
   -----------------------------------
*/
now = last = 0 ;
list[0] = u = uexp ;
mark[u] = tag ;
while ( now <= last ) {
   u = list[now++] ;
   Graph_adjAndSize(bpg->graph, u, &usize, &uadj) ;
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n    checking out u = %d : ", u) ;
      IVfp80(msgFile, usize, uadj, 20, &ierr) ;
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < usize ; ii++ ) {
      v = uadj[ii] ;
      if ( mark[v] != tag ) {
         mark[v] = tag ;
         link[v] =  u  ; 
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, 
                    "\n       adding %d with nvexp[%d] = %d to tree",
                    v, v, nvexp[v]) ;
            fflush(msgFile) ;
         }
         if ( nvexp[v] > 0 ) {
            if ( msglvl > 1 && msgFile != NULL ) {
               fprintf(msgFile, ", exposed") ;
               fflush(msgFile) ;
            }
            return(v) ;
         } else {
            IVL_listAndSize(ivl, v, &vsize, &vind) ;
            for ( jj = 0 ; jj < vsize ; jj++ ) {
               unew = vind[jj] ;
               if ( unew == -1 ) {
                  break ;
               } else if ( mark[unew] != tag ) {
                  if ( msglvl > 1 && msgFile != NULL ) {
                     fprintf(msgFile, 
                             "\n          adding %d to tree", unew) ;
                     fflush(msgFile) ;
                  }
                  mark[unew] = tag ;
                  link[unew] =  v  ;
                  list[++last] = unew ;
               }
            }
         }
      }
   }
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n leaving nonunitFindExposedNode") ;
   fflush(msgFile) ;
}
return(-1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   find the minimum number of match weights along the path from
   exposed node uexp to exposed node vexp.

   link  -- link vector, used to define alternating path, size nX + nY
   nvexp -- nvexp[u] = # of exposed vertices in u

   return value -- minimum match weight

   created -- 95oct12, cca
   -------------------------------------------------------------------
*/
static int
nonunitFindNmatch (
   BPG   *bpg,
   int   uexp,
   int   vexp,
   int   nvexp[],
   int   link[],
   IVL   *ivl,
   int   msglvl,
   FILE   *msgFile
) {
int   count, ii, msize, nmatch, u, v ;
int   *mind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL 
   || uexp < 0 || uexp > (bpg->nX + bpg->nY)
   || vexp < 0 || vexp > (bpg->nX + bpg->nY)
   || (uexp < bpg->nX && vexp < bpg->nX)
   || (uexp >= bpg->nX && vexp >= bpg->nX)
   || nvexp == NULL || link == NULL || ivl == NULL ) {
   fprintf(stderr, 
           "\n fatal error in nonunitFindNmatch(%p,%d,%d,%p,%p,%p)"
           "\n bad input\n", bpg, uexp, vexp, nvexp, link, ivl) ;
   exit(-1) ;
}
/*
   ----------------------------------------------------------
   set nmatch = min(# of exposed vertices in uexp and vexp)
   ----------------------------------------------------------
*/
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n nvexp[%d] = %d, nvexp[%d] = %d",
           uexp, nvexp[uexp], vexp, nvexp[vexp]) ;
   fflush(msgFile) ;
}
nmatch = (nvexp[uexp] <= nvexp[vexp]) ? nvexp[uexp] : nvexp[vexp] ;
if ( nmatch == 1 ) {
   return(1) ;
}
/*
   --------------------------------------------------------------- 
   find the minimum of the weight of the matched edges in the path
   --------------------------------------------------------------- 
*/
u = link[vexp] ;
while ( u != uexp ) {
   v = link[u] ;
/*
   ---------------------------------------------
   count the number of times u is matched with v
   ---------------------------------------------
*/
   IVL_listAndSize(ivl, u, &msize, &mind) ;
   for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
      if ( mind[ii] == -1 ) {
         break ;
      } else if ( mind[ii] == v ) {
         count++ ;
      }
   }
/*
   ---------------------
   check the match count
   ---------------------
*/
   if ( nmatch > count ) {
      nmatch = count ;
      if ( nmatch == 1 ) {
         return(1) ;
      }
   }
   u = link[v] ;
}

return(nmatch) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   given an alternating path, flip the edges

   uexp   -- exposed node at start of path
   vexp   -- exposed node at end of path
   nmatch -- number of matched edges to flip
   nvexp  -- nvexp[u] = # of exposed vertices in u
   link   -- link vector, used to define alternating path, 
             size nX + nY
   ivl    -- IVL object to hold matched edges

   created -- 95oct12, cca
   -------------------------------------------------------
*/
static void
nonunitFlipEdges (
   BPG   *bpg,
   int   uexp,
   int   vexp,
   int   nmatch,
   int   nvexp[],
   int   link[],
   IVL   *ivl,
   int   msglvl,
   FILE   *msgFile
) {
int   count, ierr, ii, msize, u, v, vnext ;
int   *mind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL 
   || uexp < 0 || uexp > (bpg->nX + bpg->nY)
   || vexp < 0 || vexp > (bpg->nX + bpg->nY)
   || (uexp < bpg->nX && vexp < bpg->nX)
   || (uexp >= bpg->nX && vexp >= bpg->nX)
   || nmatch <= 0
   || nvexp == NULL || link == NULL || ivl == NULL ) {
   fprintf(stderr, 
          "\n fatal error in BPG_flipMatch(%p,%d,%d,%d,%p,%p,%p)"
          "\n bad input\n", bpg, uexp, vexp, nmatch, nvexp, link, ivl) ;
   exit(-1) ;
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n match lists before edge swaps") ;
   IVL_listAndSize(ivl, vexp, &msize, &mind) ;
   fprintf(msgFile, "\n %d's match list :", vexp) ;
   IVfp80(msgFile, msize, mind, 20, &ierr) ;
   u = link[vexp] ;
   while ( u != uexp ) {
      IVL_listAndSize(ivl, u, &msize, &mind) ;
      fprintf(msgFile, "\n %d's match list :", u) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      v = link[u] ;
      IVL_listAndSize(ivl, v, &msize, &mind) ;
      fprintf(msgFile, "\n %d's match list :", v) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      u = link[v] ;
   }
   IVL_listAndSize(ivl, uexp, &msize, &mind) ;
   fprintf(msgFile, "\n %d's match list :", uexp) ;
   IVfp80(msgFile, msize, mind, 20, &ierr) ;
   fflush(msgFile) ;
}
if ( link[vexp] == uexp ) {
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n simple case, %d <---> %d", uexp, vexp) ;
      fflush(msgFile) ;
   }
/*
   ------------------------------------
   simple case, uexp --> vexp
   1. add <vexp, nmatch> to uexp's list
   2. add <uexp, nmatch> to vexp's list
   ------------------------------------
*/
   IVL_listAndSize(ivl, uexp, &msize, &mind) ;
   for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
      if ( mind[ii] == -1 ) {
         mind[ii] = vexp ;
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n   %d's list, mind[%d] = %d", 
                    uexp, ii, mind[ii]) ;
            fflush(msgFile) ;
         }
         if ( ++count == nmatch ) {
            break ;
         }
      }
   }
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n %d's match list :", uexp) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      fflush(msgFile) ;
   }
   if ( count != nmatch ) {
      fprintf(stderr, "\n fatal error in BPG_flipMatch()"
              "\n uexp = %d, vexp = %d, count = %d, nmatch = %d"
              "\n u matches with ", uexp, vexp, count, nmatch) ;
      IVfp80(stderr, msize, mind, 12, &ierr) ;
      exit(-1) ;
   }
   IVL_listAndSize(ivl, vexp, &msize, &mind) ;
   for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
      if ( mind[ii] == -1 ) {
         mind[ii] = uexp ;
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n   %d's list, mind[%d] = %d", 
                    vexp, ii, mind[ii]) ;
            fflush(msgFile) ;
         }
         if ( ++count == nmatch ) {
            break ;
         }
      }
   }
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n %d's match list :", vexp) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      fflush(msgFile) ;
   }
   if ( count != nmatch ) {
      fprintf(stderr, "\n fatal error in BPG_flipMatch()"
              "\n uexp = %d, vexp = %d, count = %d, nmatch = %d"
              "\n v matches with ", uexp, vexp, count, nmatch) ;
      IVfp80(stderr, msize, mind, 12, &ierr) ;
      exit(-1) ;
   }
} else {
/*
   -------------------------------------------------------
   uexp <--> v <==> u <--> v <==> u ... v <==> u <--> vexp
   deal with the middle edges
   -------------------------------------------------------
*/
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n complex case : %d <---> ... <---> %d",
              uexp, vexp) ;
      fflush(msgFile) ;
   }
   vnext = vexp ;
   u = link[vexp] ;
   while ( u != uexp ) {
      v = link[u] ;
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n checking out %d <===> %d <---> %d",
                 v, u, vnext) ;
         fflush(msgFile) ;
      }
/*
      ----------------------------------------------------------
      check out u's list of matches, find <v, nmatch(u,v)>
      if ( nmatch(u,v) == nmatch ) then
         replace <v, nmatch(u,v)> with <vnext, nmatch>
      else if ( nmatch(u,v) > nmatch ) then
         replace <v, nmatch(u,v)> with <v, nmatch(u,v) - nmatch>
         add <vnext, nmatch>
      else
         error
      endif
      ----------------------------------------------------------
*/
      IVL_listAndSize(ivl, u, &msize, &mind) ;
      for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
         if ( mind[ii] == v ) {
            mind[ii] = vnext ;
            if ( ++count == nmatch ) {
               break ;
            }
         }
      }
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n %d's match list :", u) ;
         IVfp80(msgFile, msize, mind, 20, &ierr) ;
         fflush(msgFile) ;
      }
      if ( count != nmatch ) {
         fprintf(stderr, "\n fatal error in BPG_flipMatch()"
                 "\n u = %d, count = %d, nmatch = %d"
                 "\n %d match list : ", u, count, nmatch, u) ;
         IVfp80(stderr, msize, mind, 12, &ierr) ;
         exit(-1) ;
      }
      IVL_listAndSize(ivl, v, &msize, &mind) ;
      for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
         if ( mind[ii] == u ) {
            mind[ii] = -1 ;
            if ( ++count == nmatch ) {
               break ;
            }
         }
      }
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n %d's match list :", v) ;
         IVfp80(msgFile, msize, mind, 20, &ierr) ;
         fflush(msgFile) ;
      }
      if ( count != nmatch ) {
         fprintf(stderr, "\n fatal error in BPG_flipMatch()"
                 "\n v = %d, count = %d, nmatch = %d"
                 "\n %d match list : ", v, count, nmatch, v) ;
         IVfp80(stderr, msize, mind, 12, &ierr) ;
         exit(-1) ;
      }
      IVL_listAndSize(ivl, vnext, &msize, &mind) ;
      for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
         if ( mind[ii] == -1 ) {
            mind[ii] = u ;
            if ( ++count == nmatch ) {
               break ;
            }
         }
      }
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n %d's match list :", vnext) ;
         IVfp80(msgFile, msize, mind, 20, &ierr) ;
         fflush(msgFile) ;
      }
      if ( count != nmatch ) {
         fprintf(stderr, "\n fatal error in BPG_flipMatch()"
                 "\n vnext = %d, count = %d, nmatch = %d"
                 "\n %d match list : ", vnext, count, nmatch, vnext) ;
         IVfp80(stderr, msize, mind, 12, &ierr) ;
         exit(-1) ;
      }
/*
      ---------------------------
      on to the next matched pair
      ---------------------------
*/
      vnext = v ;
      u = link[v] ;
   }
   IVL_listAndSize(ivl, uexp, &msize, &mind) ;
   for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
      if ( mind[ii] == -1 ) {
         mind[ii] = vnext ;
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n   %d's list, mind[%d] = %d", 
                    uexp, ii, mind[ii]) ;
            fflush(msgFile) ;
         }
         if ( ++count == nmatch ) {
            break ;
         }
      }
   }
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n %d's match list :", uexp) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      fflush(msgFile) ;
   }
   if ( count != nmatch ) {
      fprintf(stderr, "\n fatal error in BPG_flipMatch()"
              "\n uexp = %d, vnext = %d, count = %d, nmatch = %d"
              "\n u matches with ", uexp, vnext, count, nmatch) ;
      IVfp80(stderr, msize, mind, 12, &ierr) ;
      exit(-1) ;
   }
   IVL_listAndSize(ivl, vnext, &msize, &mind) ;
   for ( ii = 0, count = 0 ; ii < msize ; ii++ ) {
      if ( mind[ii] == -1 ) {
         mind[ii] = uexp ;
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n   %d's list, mind[%d] = %d", 
                    vnext, ii, mind[ii]) ;
            fflush(msgFile) ;
         }
         if ( ++count == nmatch ) {
            break ;
         }
      }
   }
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n %d's match list :", vnext) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      fflush(msgFile) ;
   }
   if ( count != nmatch ) {
      fprintf(stderr, "\n fatal error in BPG_flipMatch()"
              "\n vnext = %d, uexp = %d, count = %d, nmatch = %d"
              "\n v matches with ", vnext, uexp, count, nmatch) ;
      IVfp80(stderr, msize, mind, 12, &ierr) ;
      exit(-1) ;
   }
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n match lists after edge swaps") ;
   IVL_listAndSize(ivl, vexp, &msize, &mind) ;
   fprintf(msgFile, "\n %d's match list :", vexp) ;
   IVfp80(msgFile, msize, mind, 20, &ierr) ;
   u = link[vexp] ;
   while ( u != uexp ) {
      IVL_listAndSize(ivl, u, &msize, &mind) ;
      fprintf(msgFile, "\n %d's match list :", u) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      v = link[u] ;
      IVL_listAndSize(ivl, v, &msize, &mind) ;
      fprintf(msgFile, "\n %d's match list :", v) ;
      IVfp80(msgFile, msize, mind, 20, &ierr) ;
      u = link[v] ;
   }
   IVL_listAndSize(ivl, uexp, &msize, &mind) ;
   fprintf(msgFile, "\n %d's match list :", uexp) ;
   IVfp80(msgFile, msize, mind, 20, &ierr) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   decrement the nvexp[] values for uexp and vexp
   ----------------------------------------------
*/
nvexp[uexp] -= nmatch ;
nvexp[vexp] -= nmatch ;
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n    nvexp[%d] = %d", uexp, nvexp[uexp]) ;
   fprintf(msgFile, ", nvexp[%d] = %d", vexp, nvexp[vexp]) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   given the maximum matching, set the flags for the DM decomposition

   list    -- list vector, size max(nX, nY)
   nvexp   -- nvexp[u] = # of exposed vertices in u
   mark    -- mark vector, used to keep track of nodes in the tree
   ivl     -- IVL object to hold matched edges
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 95oct14, cca
   ------------------------------------------------------------------
*/
static void
nonunitSetFlags (
   BPG   *bpg,
   int   list[],
   int   nvexp[],
   int   mark[],
   IVL   *ivl,
   int   dmflags[], 
   int   stats[],
   int   msglvl,
   FILE   *msgFile
) {
int   ii, jj, last, now, nX, nY, x, xnew, xsize, xwght,
      y, ynew, ysize, ywght ;
int   *vwghts, *xadj, *yadj ;
/*
   ---------------
   check the input
   ---------------
*/
if (  bpg  == NULL || list == NULL || nvexp   == NULL 
   || mark == NULL || ivl  == NULL || dmflags == NULL ) {
   fprintf(stderr, 
          "\n fatal error in BPG_seDMflags(%p,%p,%p,%p,%p,%p,%p)"
          "\n bad input\n", 
          bpg, list, nvexp, mark, ivl, dmflags, stats) ;
   exit(-1) ;
}
nX     = bpg->nX ;
nY     = bpg->nY ;
vwghts = bpg->graph->vwghts ;
/*
   ----------------------------------------------
   zero the dmflags[] vector and set mark[] to -1
   ----------------------------------------------
*/
IVzero(nX + nY, dmflags) ;
IVfill(nX + nY, mark, -1) ;
IVzero(6, stats) ;
/*
   -------------------------------------
   load the list with exposed nodes in X
   -------------------------------------
*/
xwght = 0 ;
last  = -1 ;
for ( x = 0 ; x < nX ; x++ ) {
   if ( nvexp[x] > 0 ) {
      list[++last] = x ;
      mark[x]      = 1 ;
   }
   xwght += vwghts[x] ;
}
/*
   ------------------------------------------------------
   drop an alternating level structure from the exposed 
   nodes in X, and set dmflags[] for nodes in X_I and Y_E
   ------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   x = list[now++] ;
   dmflags[x] = 1 ;
   stats[0] += vwghts[x] ;
   Graph_adjAndSize(bpg->graph, x, &xsize, &xadj) ;
   for ( ii = 0 ; ii < xsize ; ii++ ) {
      y = xadj[ii] ;
      if ( mark[y] != 1 ) {
         mark[y]    = 1 ;
         dmflags[y] = 2 ;
         stats[4] += vwghts[y] ;
         IVL_listAndSize(ivl, y, &ysize, &yadj) ;
         for ( jj = 0 ; jj < ysize ; jj++ ) {
            if ( (xnew = yadj[jj]) == -1 ) {
               break ;
            } else if ( mark[xnew] != 1 ) {
               mark[xnew]   =   1  ;
               list[++last] = xnew ;
            }
         }
      }
   }
}
/*
   -------------------------------------
   load the list with exposed nodes in Y
   -------------------------------------
*/
ywght =  0 ;
last  = -1 ;
for ( y = nX ; y < nX + nY ; y++ ) {
   if ( nvexp[y] > 0 ) {
      list[++last] = y ;
      mark[y]      = 2 ;
   }
   ywght += vwghts[y] ;
}
/*
   ------------------------------------------------------
   drop an alternating level structure from the exposed 
   nodes in Y, and set dmflags[] for nodes in Y_I and X_E
   ------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   y = list[now++] ;
   dmflags[y] = 1 ;
   stats[3] += vwghts[y] ;
   Graph_adjAndSize(bpg->graph, y, &ysize, &yadj) ;
   for ( ii = 0 ; ii < ysize ; ii++ ) {
      x = yadj[ii] ;
      if ( mark[x] != 2 ) {
         mark[x]    = 2 ;
         dmflags[x] = 2 ;
         stats[1] += vwghts[x] ;
         IVL_listAndSize(ivl, x, &xsize, &xadj) ;
         for ( jj = 0 ; jj < xsize ; jj++ ) {
            if ( (ynew = xadj[jj]) == -1 ) {
               break ;
            } else if ( mark[ynew] != 2 ) {
               mark[ynew]   =   2  ;
               list[++last] = ynew ;
            }
         }
      }
   }
}
/*
   ------------------------------
   set the weights of X_R and Y_R
   ------------------------------
*/
stats[2] = xwght - stats[0] - stats[1] ;
stats[5] = ywght - stats[3] - stats[4] ;

return ; }
     
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   compute the generalized Dulmadge-Mendolsohn decomposition
   for a unit weight graph.

   bpg     -- BPG bipartite graph object
   dmflags -- flags vector
                   / 0 if x in X_R
      dmflags[x] = | 1 if x in X_I
                   \ 2 if x in X_E
                   / 0 if y in Y_R
      dmflags[y] = | 1 if y in Y_I
                   \ 2 if y in Y_E
   stats -- statistics vector
      stats[0] -- weight of X_I
      stats[1] -- weight of X_E
      stats[2] -- weight of X_R
      stats[3] -- weight of Y_I
      stats[4] -- weight of Y_E
      stats[5] -- weight of Y_R

   created -- 95oct12, cca
   ---------------------------------------------------------
*/
static void
unitDM (
   BPG   *bpg,
   int   dmflags[],
   int   stats[],
   int   msglvl,
   FILE   *msgFile
) {
int   nX, nY ;
int   *mate ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || dmflags == NULL || stats == NULL ) {
   fprintf(stderr, "\n fatal error in unitDM(%p,%p,%p)"
           "\n bad input\n", bpg, dmflags, stats) ;
   exit(-1) ;
}
nX = bpg->nX ;
nY = bpg->nY ;
mate = IVinit(nX + nY, -1) ;
if ( msglvl > 1 && msgFile != NULL ) {
   BPG_writeForHumanEye(bpg, msgFile) ;
}
/*
   -----------------------
   find a maximum matching
   -----------------------
*/
unitFindMaxMatch(bpg, mate, msglvl, msgFile) ; 
/*
   --------------------------------------
   fill the dmflags[] and stats[] vectors
   --------------------------------------
*/
unitSetFlags(bpg, mate, dmflags, stats, msglvl, msgFile) ; 

return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   find a maximum matching for a unit weight graph

   created -- 95oct14, cca
   -----------------------------------------------
*/
static void
unitFindMaxMatch (
   BPG   *bpg,
   int   mate[],
   int   msglvl,
   FILE   *msgFile
) {
int   nX, nY, rc, x, y ;
int   *link, *list, *mark ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || mate == NULL ) {
   fprintf(stderr, "\n fatal error in unitFindMaxMatch(%p,%p)"
           "\n bad input\n", bpg, mate) ;
   exit(-1) ;
}
nX = bpg->nX ;
nY = bpg->nY ;
/*
   --------------------------
   create the working storage
   --------------------------
*/
mark = IVinit(nX + nY, -1) ;
link = IVinit(nX + nY, -1) ;
list = IVinit(nX + nY, -1) ;
/*
   ---------------------------------------------------------
   look for an augmenting path rooted at each exposed vertex
   ---------------------------------------------------------
*/
if ( nX <= nY ) {
   for ( x = 0 ; x < nX ; x++ ) {
      if ( mate[x] == -1 ) {
         rc = unitAugmentingPath(bpg, x, mate, mark, link, list,
                                 msglvl, msgFile) ; 
      }
   }
} else {
   for ( y = nX ; y < nX + nY ; y++ ) {
      if ( mate[y] == -1 ) {
         rc = unitAugmentingPath(bpg, y, mate, mark, link, list,
                                 msglvl, msgFile) ; 
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(link) ;
IVfree(list) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   look for a augmenting path starting at node uexp.
   if one found, flip the match edges.

   return value -- 1 if success
                   0 otherwise

   created -- 95oct14, cca
   -------------------------------------------------
*/
static int
unitAugmentingPath (
   BPG   *bpg,
   int   uexp,
   int   mate[],
   int   mark[],
   int   link[], 
   int   list[],
   int   msglvl,
   FILE   *msgFile
) {
int   ii, last, nextv, now, u, usize, v ;
int   *uadj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( bpg == NULL || uexp < 0 || (bpg->nX + bpg->nY <= uexp) 
   || mate == NULL || mark == NULL || link == NULL || list == NULL 
   || mate[uexp] != -1 ) {
   fprintf(stderr, 
           "\n fatal error in unitAugmentingPath(%p,%d,%p,%p,%p,%p)"
           "\n bad input\n", bpg, uexp, mate, mark, link, list) ;
   exit(-1) ;
}
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n\n uexp = %d", uexp) ;
}
now = last = 0 ;
list[0]    = uexp ;
mark[uexp] = uexp ;
while ( now <= last ) {
   u = list[now++] ;
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n    checking out %d", u) ;
   }
   Graph_adjAndSize(bpg->graph, u, &usize, &uadj) ;
   for ( ii = 0 ; ii < usize ; ii++ ) {
      v = uadj[ii] ;
      if ( mark[v] != uexp ) {
/*
         ------------------------------------
         v is not yet in the alternating tree
         ------------------------------------
*/
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n       adding v = %d to tree", v) ;
         }
         mark[v] = uexp ;
         link[v] =   u  ;
         if ( mate[v] == -1 ) {
            if ( msglvl > 1 && msgFile != NULL ) {
               fprintf(msgFile, ", exposed") ;
            }
/*
            ---------------------------------------------
            v is exposed, flip the match edges and return
            ---------------------------------------------
*/
            while ( v != -1 ) {
               u = link[v] ;
               nextv = mate[u] ;
               if ( msglvl > 1 && msgFile != NULL ) {
                  fprintf(msgFile, "\n       setting %d <===> %d", v, u);
               }
               mate[u] = v ;
               mate[v] = u ;
               v = nextv ;
            }
            return(1) ;
         } else {
/*
            ------------------------
            put v's mate on the list
            ------------------------
*/
            if ( msglvl > 1 && msgFile != NULL ) {
               fprintf(msgFile, 
                       "\n       adding u = %d to tree", mate[v]) ;
            }
            list[++last] = mate[v] ;
         }
      }
   }
}
return(0) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   set the dmflags[] and stats[] arrays for a unit weight graph
   ------------------------------------------------------------
*/
static void
unitSetFlags (
   BPG   *bpg,
   int   mate[],
   int   dmflags[], 
   int   stats[],
   int   msglvl,
   FILE   *msgFile
) {
int   ierr, ii, last, now, nX, nY, x, xmate, xsize, y, ymate, ysize ;
int   *list, *xadj, *yadj ;
/*
   ---------------
   check the input
   ---------------
*/
if (  bpg == NULL || mate == NULL || dmflags == NULL || stats == NULL ){
   fprintf(stderr, "\n fatal error in unitSetFlags(%p,%p,%p,%p)"
           "\n bad input\n", bpg, mate, dmflags, stats) ;
   exit(-1) ; 
}
nX = bpg->nX ;
nY = bpg->nY ;
list = IVinit(nX + nY, -1) ;
/*
   --------------------------------------
   zero the dmflags[] and stats[] vectors
   --------------------------------------
*/
IVzero(nX + nY, dmflags) ;
IVzero(6, stats) ;
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n MATE") ;
   IVfp80(msgFile, nX + nY, mate, 80, &ierr) ;
}
/*
   -------------------------------------
   load the list with exposed nodes in X
   -------------------------------------
*/
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n\n exposed nodes in X") ;
}
last  = -1 ;
for ( x = 0 ; x < nX ; x++ ) {
   if ( mate[x] == -1 ) {
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n loading x = %d", x) ;
      }
      list[++last] = x ;
      dmflags[x]   = 1 ;
   }
}
/*
   ------------------------------------------------------
   drop an alternating level structure from the exposed 
   nodes in X, and set dmflags[] for nodes in X_I and Y_E
   ------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   x = list[now++] ;
   Graph_adjAndSize(bpg->graph, x, &xsize, &xadj) ;
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n adj(%d) :", x) ;
      IVfp80(msgFile, xsize, xadj, 10, &ierr) ;
      fflush(msgFile) ;
   }
   for ( ii = 0 ; ii < xsize ; ii++ ) {
      y = xadj[ii] ;
      if ( dmflags[y] == 0 ) {
         if ( msglvl > 1 && msgFile != NULL ) {
            fprintf(msgFile, "\n    adding y = %d", y) ;
         }
         dmflags[y] = 2 ;
         if ( (xmate = mate[y]) == -1 ) {
            fprintf(stderr, "\n fatal error in unitSetFlags"
                    "\n now = %d, x = %d, y = %d, mate = -1", 
                    now, x, y) ;
            fprintf(stderr, "\n adj(%d) :", x) ;
            IVfp80(stderr, xsize, xadj, 10, &ierr) ;
            exit(-1) ;
         } else if ( dmflags[xmate] == 0 ) {
            if ( msglvl > 1 && msgFile != NULL ) {
               fprintf(msgFile, "\n    adding xmate = %d", xmate) ;
            }
            dmflags[xmate] =   1   ;
            list[++last]   = xmate ;
         }
      }
   }
}
/*
   -------------------------------------
   load the list with exposed nodes in Y
   -------------------------------------
*/
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, "\n\n exposed nodes in X") ;
}
last  = -1 ;
for ( y = 0 ; y < nX + nY ; y++ ) {
   if ( mate[y] == -1 ) {
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n loading y = %d", y) ;
      }
      list[++last] = y ;
      dmflags[y]   = 1 ;
   }
}
/*
   ------------------------------------------------------
   drop an alternating level structure from the exposed 
   nodes in X, and set dmflags[] for nodes in X_I and Y_E
   ------------------------------------------------------
*/
now = 0 ;
while ( now <= last ) {
   y = list[now++] ;
   Graph_adjAndSize(bpg->graph, y, &ysize, &yadj) ;
   for ( ii = 0 ; ii < ysize ; ii++ ) {
      x = yadj[ii] ;
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n    loading x = %d", x) ;
      }
      if ( dmflags[x] == 0 ) {
         dmflags[x] = 2 ;
         if ( (ymate = mate[x]) == -1 ) {
            fprintf(stderr, "\n fatal error in unitSetFlags"
                    "\n x = %d, mate = -1", x) ;
            exit(-1) ;
         } else if ( dmflags[ymate] == 0 ) {
            if ( msglvl > 1 && msgFile != NULL ) {
               fprintf(msgFile, "\n    loading ymate = %d", ymate) ;
            }
            dmflags[ymate] =   1   ;
            list[++last]   = ymate ;
         }
      }
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(list) ;
/*
   ------------------------------
   set the weights of X_R and Y_R
   ------------------------------
*/
stats[0] = stats[1] = stats[2] = stats[3] = stats[4] = stats[5] = 0 ;
for ( x = 0 ; x < nX ; x++ ) {
   switch ( dmflags[x] ) {
   case 0 : stats[2]++ ; break ;
   case 1 : stats[0]++ ; break ;
   case 2 : stats[1]++ ; break ;
   }
}
for ( y = nX ; y < nX + nY ; y++ ) {
   switch ( dmflags[y] ) {
   case 0 : stats[5]++ ; break ;
   case 1 : stats[3]++ ; break ;
   case 2 : stats[4]++ ; break ;
   }
}

return ; }
     
/*--------------------------------------------------------------------*/
