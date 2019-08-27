/*  identifyWideSep.c  */

#include "../GPart.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------
   identify the wide separator
 
   return -- IV object that holds the nodes in the wide separator

   created -- 96oct21, cca
   --------------------------------------------------------------
*/
IV *
GPart_identifyWideSep (
   GPart   *gpart,
   int     nlevel1,
   int     nlevel2
) {
FILE    *msgFile ;
Graph   *g ;
int     count, first, ierr, ii, ilevel, last, msglvl,
        nfirst, now, nsecond, nsep, nvtx, v, vsize, w ;
int     *compids, *list, *mark, *vadj ;
IV      *sepIV ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL || (g = gpart->g) == NULL 
   || nlevel1 < 0 || nlevel2 < 0 ) {
  fprintf(stderr, "\n fatal error in GPart_identifyWideSep(%p,%d,%d)"
           "\n bad input\n", gpart, nlevel1, nlevel2) ;
   exit(-1) ;
}
g       = gpart->g ;
compids = IV_entries(&gpart->compidsIV) ;
nvtx    = g->nvtx ;
mark    = IVinit(nvtx, -1) ;
list    = IVinit(nvtx, -1) ;
msglvl  = gpart->msglvl ;
msgFile = gpart->msgFile ;
/*
   --------------------------------------
   load the separator nodes into the list
   --------------------------------------
*/
nsep = 0 ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == 0 ) {
      list[nsep++] = v ;
      mark[v] = 0 ;
   }
}
count = nsep ;
if ( msglvl > 1 ) {
   fprintf(msgFile, 
           "\n GPart_identifyWideSep : %d separator nodes loaded", 
           count) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   IVfp80(msgFile, nsep, list, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   loop over the number of levels out that form 
   the wide separator towards the first component
   ----------------------------------------------
*/
if ( nlevel1 >= 1 ) {
   first = count ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n level = %d, first = %d", 1, first) ;
      fflush(msgFile) ;
   }
   for ( now = 0 ; now < nsep ; now++ ) {
      v = list[now] ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n %d : ", v) ;
         IVfp80(msgFile, vsize, vadj, 80, &ierr) ;
         fflush(msgFile) ;
      }
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < nvtx && mark[w] == -1 && compids[w] == 1 ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n    adding %d to list", w) ;
               fflush(msgFile) ;
            }
            list[count++] = w ;
            mark[w] = 1 ;
         }
      }
   }
   now = first ;
   for ( ilevel = 2 ; ilevel <= nlevel1 ; ilevel++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n level = %d, first = %d", ilevel, first);
         fflush(msgFile) ;
      }
      last = count - 1 ;
      while ( now <= last ) {
         v = list[now++] ;
         Graph_adjAndSize(g, v, &vsize, &vadj) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n %d : ", v) ;
            IVfp80(msgFile, vsize, vadj, 80, &ierr) ;
            fflush(msgFile) ;
         }
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( w < nvtx && mark[w] == -1 && compids[w] == 1 ) {
               if ( msglvl > 2 ) {
                  fprintf(msgFile, "\n    adding %d to list", w) ;
                  fflush(msgFile) ;
               }
               mark[w] = 1 ;
               list[count++] = w ;
            }
         }
      }
   }
}
nfirst = count - nsep ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n %d nodes added from the first component", nfirst) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   IVfp80(msgFile, nfirst, &list[nsep], 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   loop over the number of levels out that form 
   the wide separator towards the second component
   ----------------------------------------------
*/
if ( nlevel2 >= 1 ) {
   first = count ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n level = %d, first = %d", 1, first) ;
      fflush(msgFile) ;
   }
   for ( now = 0 ; now < nsep ; now++ ) {
      v = list[now] ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n %d : ", v) ;
         IVfp80(msgFile, vsize, vadj, 80, &ierr) ;
         fflush(msgFile) ;
      }
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < nvtx && mark[w] == -1 && compids[w] == 2 ) {
            if ( msglvl > 2 ) {
               fprintf(msgFile, "\n    adding %d to list", w) ;
               fflush(msgFile) ;
            }
            list[count++] = w ;
            mark[w] = 2 ;
         }
      }
   }
   now = first ;
   for ( ilevel = 2 ; ilevel <= nlevel2 ; ilevel++ ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n\n level = %d, first = %d", ilevel, first);
         fflush(msgFile) ;
      }
      last = count - 1 ;
      while ( now <= last ) {
         v = list[now++] ;
         Graph_adjAndSize(g, v, &vsize, &vadj) ;
         if ( msglvl > 2 ) {
            fprintf(msgFile, "\n %d : ", v) ;
            IVfp80(msgFile, vsize, vadj, 80, &ierr) ;
            fflush(msgFile) ;
         }
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( w < nvtx && mark[w] == -1 && compids[w] == 2 ) {
               if ( msglvl > 2 ) {
                  fprintf(msgFile, "\n    adding %d to list", w) ;
                  fflush(msgFile) ;
               }
               mark[w] = 2 ;
               list[count++] = w ;
            }
         }
      }
   }
}
nsecond = count - nsep - nfirst ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n %d nodes added from the second component", 
           nsecond) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   IVfp80(msgFile, nsecond, &list[nsep + nfirst], 80, &ierr) ;
   fflush(msgFile) ;
}
IVqsortUp(count, list) ;
/*
   --------------------
   create the IV object
   --------------------
*/
sepIV = IV_new() ;
IV_init(sepIV, count, NULL) ;
IVcopy(count, IV_entries(sepIV), list) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n separator has %d nodes", IV_size(sepIV)) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n sepIV") ;
   IV_writeForHumanEye(sepIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(mark) ;
IVfree(list) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n return from GPart_identifyWideSep") ;
   fflush(msgFile) ;
}
 
return(sepIV) ; }

/*--------------------------------------------------------------------*/
