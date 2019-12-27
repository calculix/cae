/*  smoothYSep.c  */

#include "../GPart.h"
#include "../../Network.h"

/*--------------------------------------------------------------------*/
/*
   -----------------
   static prototypes
   -----------------
*/
static Network * 
createNetwork ( Graph *g, int compids[], IV *YVmapIV, IV *YCmapIV, 
                IV *NYmapIV, int msglvl, FILE *msgFile ) ;
static void getNewCompids ( int nnode, int NYmap[], int YCmap[],
                            int mark[], int Ycompids[],
                            int   msglvl, FILE  *msgFile ) ;
static float eval ( float alpha, float wS, float wB, float wW) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   smooth the wide separator given by YVmapIV by forming 
   the associated network, solving the max flow problem, 
   and examining two min-cuts.

   YVmapIV -- map from wide vertices Y to vertices V
   YCmapIV -- map from wide vertices Y to {0,1,2,3}
     YCmap[y] = 0 --> treat y as an internal node
                      adjacent to neither component
     YCmap[y] = 1 --> treat y as adjacent only to component 1
     YCmap[y] = 2 --> treat y as adjacent only to component 2
     YCmap[y] = 3 --> treat y as adjacent to components 1 and 2
   alpha -- cost function parameter

   created -- 96jun08, cca
   ------------------------------------------------------------
*/
float
GPart_smoothYSep (
   GPart   *gpart,
   IV      *YVmapIV,
   IV      *YCmapIV,
   float   alpha
) {
FILE      *msgFile ;
float     bestcost, newcost1, newcost2 ;
Graph     *g ;
Ideq      *deq ;
int       ierr, msglvl, nnode, newcomp, nY, oldcomp, v, vwght, y ;
int       *compids, *cweights, *mark, *vwghts, 
          *Ycompids1, *Ycompids2, *YCmap, *YVmap ;
int       deltas[3] ;
IV        *NYmapIV ;
Network   *network ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL || (g = gpart->g) == NULL 
   || YVmapIV == NULL 
   || (nY = IV_size(YVmapIV)) <= 0
   || (YVmap = IV_entries(YVmapIV)) == NULL
   || YCmapIV == NULL 
   || IV_size(YCmapIV) != nY
   || (YCmap = IV_entries(YCmapIV)) == NULL
   || alpha < 0.0 ) {
   fprintf(stderr, "\n fatal error in GPart_smoothSep(%p,%p,%p,%f)"
           "\n bad input\n", gpart, YVmapIV, YCmapIV, alpha) ;
   exit(-1) ;
}
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
vwghts   = g->vwghts ;
msgFile  = gpart->msgFile ;
msglvl   = gpart->msglvl  ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n YVmapIV") ;
   IV_writeForHumanEye(YVmapIV, msgFile) ;
   fprintf(msgFile, "\n YCmapIV") ;
   IV_writeForHumanEye(YCmapIV, msgFile) ;
}
/*
   ------------------
   create the network
   ------------------
*/
NYmapIV = IV_new() ;
network = createNetwork(g, compids, YVmapIV, YCmapIV, NYmapIV,
                        msglvl, msgFile) ;
Network_setMessageInfo(network, msglvl, msgFile) ;
nnode   = network->nnode ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n network : %d nodes, %d arcs", 
           nnode, network->narc) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   Network_writeForHumanEye(network, msgFile) ;
}
/*
   --------------------
   solve for a max flow
   --------------------
*/
Network_findMaxFlow(network) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n after max flow solution, %d arc traversals",
           network->ntrav) ;
   fflush(msgFile) ;
}
if ( msglvl > 4 ) {
   Network_writeForHumanEye(network, msgFile) ;
}
/*
   ------------------------------
   evaluate the present partition
   ------------------------------
*/
bestcost = eval(alpha, cweights[0], cweights[1], cweights[2]) ;
if ( msglvl > 1 ) {
   int   maxval, minval, val1, val2 ;
   fprintf(msgFile, "\n present partition : < %d, %d, %d >",
           cweights[0], cweights[1], cweights[2]) ;
   val1 =  cweights[1] ;
   val2 =  cweights[2] ;
   minval = (val1 <= val2) ? val1 : val2 ;
   maxval = (val1 > val2) ? val1 : val2 ;
   if ( minval == 0 ) {
      fprintf(msgFile, ", imbalance = infinite") ;
   } else {
      fprintf(msgFile, 
              ", imbalance = %6.3f", ((double) maxval)/minval) ;
   }
   fprintf(msgFile, ", cost = %.2f", bestcost) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   set up some working storage
   ---------------------------
*/
deq = Ideq_new() ;
Ideq_resize(deq, nnode) ;
mark      = IVinit(nnode, -1) ;
Ycompids1 = IVinit(nY, -1) ;
Ycompids2 = IVinit(nY, -1) ;
/*
   ---------------------------------------
   find a min-cut starting from the source
   and evaluate the partition
   ---------------------------------------
*/
Network_findMincutFromSource(network, deq, mark) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n mark starting from source") ;
   IVfp80(msgFile, nnode, mark, 80, &ierr) ;
   fflush(msgFile) ;
}
getNewCompids(nnode, IV_entries(NYmapIV), YCmap, mark, Ycompids1,
              msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n Ycompids1") ;
   IVfp80(msgFile, nY, Ycompids1, 80, &ierr) ;
   fflush(msgFile) ;
}
deltas[0] = deltas[1] = deltas[2] = 0 ;
for ( y = 0 ; y < nY ; y++ ) {
   v       = YVmap[y] ;
   oldcomp = compids[v] ;
   newcomp = Ycompids1[y] ;
   if ( msglvl > 5 ) {
      fprintf(msgFile, "\n y = %d, v = %d, oldcomp = %d, newcomp = %d",
              y, v, oldcomp, newcomp) ;
      fflush(msgFile) ;
   }
   if ( oldcomp != newcomp ) {
      vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
      deltas[oldcomp] -= vwght ;
      deltas[newcomp] += vwght ;
   }
}
newcost1 = eval(alpha, cweights[0] + deltas[0],
                cweights[1] + deltas[1], cweights[2] + deltas[2]) ;
if ( msglvl > 1 ) {
   int   maxval, minval, val1, val2 ;
   fprintf(msgFile, 
      "\n min-cut from source: < %d, %d, %d >",
        cweights[0] + deltas[0], cweights[1] + deltas[1], 
        cweights[2] + deltas[2]) ;
   val1 =  cweights[1] + deltas[1] ;
   val2 =  cweights[2] + deltas[2] ;
   minval = (val1 <= val2) ? val1 : val2 ;
   maxval = (val1 > val2) ? val1 : val2 ;
   if ( minval == 0 ) {
      fprintf(msgFile, ", imbalance = infinite") ;
   } else {
      fprintf(msgFile, 
              ", imbalance = %6.3f", ((double) maxval)/minval) ;
   }
   fprintf(msgFile, ", newcost1 = %.2f", newcost1) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   find a min-cut starting from the sink
   and evaluate the partition
   ---------------------------------------
*/
Network_findMincutFromSink(network, deq, mark) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n mark starting from sink") ;
   IVfp80(msgFile, nnode, mark, 80, &ierr) ;
   fflush(msgFile) ;
}
getNewCompids(nnode, IV_entries(NYmapIV), YCmap, mark, Ycompids2,
              msglvl, msgFile) ;
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n Ycompids2") ;
   IVfp80(msgFile, nY, Ycompids2, 80, &ierr) ;
   fflush(msgFile) ;
}
deltas[0] = deltas[1] = deltas[2] = 0 ;
for ( y = 0 ; y < nY ; y++ ) {
   v       = YVmap[y] ;
   oldcomp = compids[v] ;
   newcomp = Ycompids2[y] ;
   if ( oldcomp != newcomp ) {
      vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
      deltas[oldcomp] -= vwght ;
      deltas[newcomp] += vwght ;
   }
}
newcost2 = eval(alpha, cweights[0] + deltas[0],
                cweights[1] + deltas[1], cweights[2] + deltas[2]) ;
if ( msglvl > 1 ) {
   int   maxval, minval, val1, val2 ;
   fprintf(msgFile, 
      "\n min-cut from sink: < %d, %d, %d >",
        cweights[0] + deltas[0], cweights[1] + deltas[1], 
        cweights[2] + deltas[2]) ;
   val1 =  cweights[1] + deltas[1] ;
   val2 =  cweights[2] + deltas[2] ;
   minval = (val1 <= val2) ? val1 : val2 ;
   maxval = (val1 > val2) ? val1 : val2 ;
   if ( minval == 0 ) {
      fprintf(msgFile, ", imbalance = infinite") ;
   } else {
      fprintf(msgFile, 
              ", imbalance = %6.3f", ((double) maxval)/minval) ;
   }
   fprintf(msgFile, ", newcost1 = %.2f", newcost2) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   decide whether to accept either new partition
   ---------------------------------------------
*/
if ( newcost1 <= newcost2 && newcost1 < bestcost ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, 
              "\n accepting new partition from source min-cut") ;
      fflush(msgFile) ;
   }
   for ( y = 0 ; y < nY ; y++ ) {
      v       = YVmap[y] ;
      oldcomp = compids[v] ;
      newcomp = Ycompids1[y] ;
      if ( oldcomp != newcomp ) {
         compids[v] = newcomp ;
         vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
         cweights[oldcomp] -= vwght ;
         cweights[newcomp] += vwght ;
      }
   }
   bestcost = newcost1 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", cost = %.3f", bestcost) ;
      fflush(msgFile) ;
   }
} else if ( newcost2 <= newcost1 && newcost2 < bestcost ) {
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n accepting new partition from sink min-cut") ;
      fflush(msgFile) ;
   }
   for ( y = 0 ; y < nY ; y++ ) {
      v       = YVmap[y] ;
      oldcomp = compids[v] ;
      newcomp = Ycompids2[y] ;
      if ( oldcomp != newcomp ) {
         compids[v] = newcomp ;
         vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
         cweights[oldcomp] -= vwght ;
         cweights[newcomp] += vwght ;
      }
   }
   bestcost = newcost2 ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, ", cost = %.3f", bestcost) ;
      fflush(msgFile) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Network_free(network) ;
IV_free(NYmapIV) ;
IVfree(mark) ;
IVfree(Ycompids1) ;
IVfree(Ycompids2) ;
Ideq_free(deq) ;

return(bestcost) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------
   create the network
   ------------------
*/
static Network *
createNetwork (
   Graph   *g,
   int     compids[],
   IV      *YVmapIV,
   IV      *YCmapIV,
   IV      *NYmapIV,
   int     msglvl,
   FILE    *msgFile
) {
int       first, ierr, ii, i0, i1, i12, i2, maxcap, mnode, 
          n0, n1, n12, n2, nvtx, nY, second, sink, source, v, vnet,
          vsize, vwght, w, wnet, y, zwid ;
int       *NYmap, *vadj, *VNmap, *vwghts, *YCmap, *YVmap ;
Network   *network ;
/*
   ---------------
   check the input
   ---------------
*/
if ( g == NULL || (nvtx = g->nvtx) <= 0 || compids == NULL
     || YVmapIV == NULL 
     || (nY = IV_size(YVmapIV)) <= 0
     || (YVmap = IV_entries(YVmapIV)) == NULL
     || YCmapIV == NULL 
     || nY != IV_size(YCmapIV)
     || (YCmap = IV_entries(YCmapIV)) == NULL
     || NYmapIV == NULL ) {
   fprintf(stderr, "\n fatal error in createNetwork(%p,%p,%p,%p,%p)"
           "\n bad input\n", g, compids, YVmapIV, NYmapIV, YCmapIV) ;
   exit(-1) ;
}
vwghts = g->vwghts ;
if ( vwghts == NULL ) {
   maxcap = nvtx ;
} else {
   maxcap = IVsum(nvtx, vwghts) ;
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n maxcap = %d", maxcap) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------
   count the number of nodes in each of the four sets
   --------------------------------------------------
*/
n1 = n12 = n0 = n2 = 0 ;
for ( y = 0 ; y < nY ; y++ ) {
   switch ( YCmap[y] ) {
   case 0  : n0++  ; break ;
   case 1  : n1++  ; break ;
   case 2  : n2++  ; break ;
   case 3  : n12++ ; break ;
   default : 
      fprintf(stderr, "\n fatal error, y = %d, YCmap[%d] = %d",
              y, y, YCmap[y]) ;
      exit(-1) ;
   }
}
mnode = 1 + n1 + n12 + 2*n0 + n2 + 1 ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n n12 = %d, n1 = %d, n0 = %d, n2 = %d", 
           n12, n1, n0, n2) ;
   fprintf(msgFile, "\n %d nodes in the network", mnode) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------
   set up the N --> Y and V --> N maps
   -----------------------------------
*/
IV_init(NYmapIV, mnode, NULL) ;
NYmap = IV_entries(NYmapIV) ;
VNmap = IVinit(nvtx, -1) ;
i12 = 1 ;
i1  = 1 + n12 ;
i0  = 1 + n12 + n1 ;
i2  = 1 + n12 + n1 + 2*n0 ;
for ( y = 0 ; y < nY ; y++ ) {
   v = YVmap[y] ;
   switch ( YCmap[y] ) {
   case 0 :
      NYmap[i0]   =  y ;
      NYmap[i0+1] =  y ;
      VNmap[v]    = i0 ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n comp 0 : y = %d, v = %d, i0 = %d and %d", 
                 y, v, i0, i0+1) ;
         fflush(msgFile) ;
      }
      i0 += 2 ;
      break ;
   case 1 :
      NYmap[i1] =  y ;
      VNmap[v]  = i1 ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n comp 1 : y = %d, v = %d, i1 = %d", 
                 y, v, i1) ;
         fflush(msgFile) ;
      }
      i1++ ;
      break ;
   case 2 :
      NYmap[i2] =  y ;
      VNmap[v]  = i2 ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n comp 2 : y = %d, v = %d, i2 = %d", 
                 y, v, i2) ;
         fflush(msgFile) ;
      }
      i2++ ;
      break ;
   case 3 :
      NYmap[i12] =  y  ;
      VNmap[v]   = i12 ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n comp 12 : y = %d, v = %d, i12 = %d", 
                 y, v, i12) ;
         fflush(msgFile) ;
      }
      i12++ ;
      break ;
   default :
      fprintf(stderr,
              "\n fatal error, y = %d, v = %d, YCmap[%d] = %d",
              y, v, y, YCmap[y]) ;
      exit(-1) ;
   }
}
if ( msglvl > 4 ) {
   fprintf(msgFile, "\n NYmapIV") ;
   IV_writeForHumanEye(NYmapIV, msgFile) ;
   fprintf(msgFile, "\n VNmap") ;
   IVfp80(msgFile, nvtx, VNmap, 80, &ierr) ;
         fflush(msgFile) ;
}
/*
   -------------------------
   create the Network object
   -------------------------
*/
source  = 0 ;
sink    = mnode - 1 ;
network = Network_new() ;
Network_init(network, mnode, 0) ;
/*
   ---------------
   insert the arcs
   ---------------
*/
for ( y = 0 ; y < nY ; y++ ) {
   v    = YVmap[y] ;
   vnet = VNmap[v] ;
   vwght  = (vwghts != NULL) ? vwghts[v] : 1 ;
   if ( msglvl > 4 ) {
      fprintf(msgFile, 
              "\n checking out y = %d, v = %d, vnet = %d, vwght = %d",
              y, v, vnet, vwght) ;
      fflush(msgFile) ;
   }
   switch ( YCmap[y] ) {
   case 0 :
      first  = vnet ;
      second = first + 1 ;
      Network_addArc(network, first, second, vwght, 0) ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n S: arc (%d, %d) --> (%d,%d)",
                 v, v, first, second) ;
         fflush(msgFile) ;
      }
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < nvtx && v != w && (wnet = VNmap[w]) != -1 ) {
            zwid = NYmap[wnet] ;
            if ( YCmap[zwid] == 0 ) {
               first  = vnet + 1 ;
               second = wnet ;
               if ( msglvl > 4 ) {
                  fprintf(msgFile, "\n S: arc (%d, %d) --> (%d,%d)",
                          v, w, first, second) ;
                  fflush(msgFile) ;
               }
               Network_addArc(network, first, second, maxcap, 0) ;
            }
         }
      }
      break ;
   case 1 :
      first  = source ;
      second = vnet ;
      Network_addArc(network, first, second, vwght, 0) ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n B: arc (source, %d) --> (%d,%d)",
                 v, first, second) ;
         fflush(msgFile) ;
      }
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( msglvl > 4 ) {
            fprintf(msgFile, "\n w = %d", w) ;
            fflush(msgFile) ;
         }
         if ( w < nvtx && v != w && (wnet = VNmap[w]) != -1 ) {
            zwid = NYmap[wnet] ;
            if ( msglvl > 4 ) {
               fprintf(msgFile, "\n wnet = %d", wnet) ;
               fflush(msgFile) ;
               fprintf(msgFile, "\n zwid = %d", zwid) ;
               fflush(msgFile) ;
            }
            if ( YCmap[zwid] != 1 ) {
               first  = vnet ;
               second = wnet ;
               if ( msglvl > 4 ) {
                  fprintf(msgFile, "\n B: arc (%d, %d) --> (%d,%d)",
                       v, w, first, second) ;
                  fflush(msgFile) ;
               }
               Network_addArc(network, first, second, maxcap, 0) ;
            }
         }
      }
      break ;
   case 2 :
      first  = vnet ;
      second = sink ;
      Network_addArc(network, first, second, vwght, 0) ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n B: arc (%d, sink) --> (%d,%d)",
                 v, first, second) ;
         fflush(msgFile) ;
      }
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < nvtx && v != w && (wnet = VNmap[w]) != -1 ) {
            zwid = NYmap[wnet] ;
            if ( YCmap[zwid] == 0 ) {
               first  = wnet + 1 ;
               second = vnet ;
               if ( msglvl > 4 ) {
                  fprintf(msgFile, "\n B: arc (%d, %d) --> (%d,%d)",
                          w, v, first, second) ;
                  fflush(msgFile) ;
               }
               Network_addArc(network, first, second, maxcap, 0) ;
            }
         }
      }
      break ;
   case 3 :
      first  = source ;
      second = vnet ;
      Network_addArc(network, first, second, vwght, 0) ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n B: arc (source, %d) --> (%d,%d)",
                 v, first, second) ;
         fflush(msgFile) ;
      }
      first  = vnet ;
      second = sink ;
      Network_addArc(network, first, second, vwght, 0) ;
      if ( msglvl > 4 ) {
         fprintf(msgFile, "\n B: arc (%d, sink) --> (%d,%d)",
                 v, first, second) ;
         fflush(msgFile) ;
      }
      break ;
   default :
      fprintf(stderr,
              "\n fatal error, y = %d, v = %d, YCmap[%d] = %d",
              y, v, y, YCmap[y]) ;
      exit(-1) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(VNmap) ;

return(network) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------
   partition evaluation function
   -----------------------------
*/
static float
eval (
   float   alpha,
   float   wS,
   float   wB,
   float   wW
) {
float   cost ;
if ( wB == 0 || wW == 0 ) {
   cost = (wS + wB + wW) * (wS + wB + wW) ;
} else if ( wB >= wW ) {
   cost = wS*(1. + (alpha*wB)/wW) ;
} else {
   cost = wS*(1. + (alpha*wW)/wB) ;
}
return(cost) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------
   given a min-cut via the mark[] vector,
   get the new component ids

   created -- 96jun08, cca
   --------------------------------------
*/
static void
getNewCompids (
   int   nnode,
   int   NYmap[],
   int   YCmap[],
   int   mark[],
   int   Ycompids[],
   int   msglvl,
   FILE  *msgFile
) {
int   sink, y, ynet ;
/*
   ------------------------------------------
   decide whether to accept the new separator
   ------------------------------------------
*/
sink = nnode - 1 ;
ynet = 1 ;
while ( ynet < sink ) {
   y = NYmap[ynet] ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n ynet = %d, y = %d, YCmap[%d] = %d",
              ynet, y, y, YCmap[y]) ;
      fflush(msgFile) ;
   }
   switch ( YCmap[y] ) {
   case 0 :
      if ( mark[ynet] != mark[ynet+1] ) {
         Ycompids[y] = 0 ;
      } else { 
         Ycompids[y] = mark[ynet] ;
      }
      ynet += 2 ;
      break ;
   case 1 :
      if ( mark[ynet] == 1 ) {
         Ycompids[y] = 1 ;
      } else {
         Ycompids[y] = 0 ;
      }
      ynet++ ;
      break ;
   case 2 :
      if ( mark[ynet] == 2 ) {
         Ycompids[y] = 2 ;
      } else {
         Ycompids[y] = 0 ;
      }
      ynet++ ;
      break ;
   case 3 :
      Ycompids[y] = 0 ;
      ynet++ ;
      break ;
   default :
      fprintf(stderr, "\n fatal error in getNewCompids()"
              "\n ynet = %d, y = %d, YCmap[%d] = %d",
              ynet, y, y, YCmap[y]) ;
      exit(-1) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, ", Ycompids[%d] = %d", y, Ycompids[y]) ;
      fflush(msgFile) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
