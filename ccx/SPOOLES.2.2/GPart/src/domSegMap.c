/*  domSegMap.c  */

#include "../GPart.h"
#include "../../timings.h"
#define  NTIMES 12
static int      icputimes ;
static double   cputimes[NTIMES] ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------------------
   fill *pndom with ndom, the number of domains.
   fill *pnseg with nseg, the number of segments.
   domains are numbered in [0, ndom), segments in [ndom,ndom+nseg).

   return -- an IV object that contains the map 
             from vertices to segments

   created -- 99feb29, cca
   --------------------------------------------------------------------
*/
IV *
GPart_domSegMap (
   GPart   *gpart,
   int     *pndom,
   int     *pnseg
) {
FILE    *msgFile ;
Graph   *g ;
int     adjdom, count, d, first, ierr, ii, jj1, jj2, last, ndom, 
        msglvl, nextphi, nPhi, nPsi, nV, phi, phi0, phi1, phi2, phi3, 
        psi, sigma, size, size0, size1, size2, v, vsize, w ;
int     *adj, *adj0, *adj1, *adj2, *compids, *dmark, *dsmap, *head, 
        *link, *list, *offsets, *PhiToPsi, *PhiToV, *PsiToSigma, 
        *vadj, *VtoPhi ;
IV      *dsmapIV ;
IVL     *PhiByPhi, *PhiByPowD, *PsiByPowD ;
/*
   --------------------
   set the initial time
   --------------------
*/
icputimes = 0 ;
MARKTIME(cputimes[icputimes]) ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL 
   || (g = gpart->g) == NULL 
   || pndom == NULL
   || pnseg == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_domSegMap(%p,%p,%p)"
           "\n bad input\n", gpart, pndom, pnseg) ;
   exit(-1) ;
}
compids = IV_entries(&gpart->compidsIV) ;
msglvl  = gpart->msglvl  ;
msgFile = gpart->msgFile ;
/*
   ------------------------
   create the map IV object
   ------------------------
*/
nV = g->nvtx ;
dsmapIV = IV_new() ;
IV_init(dsmapIV, nV, NULL) ;
dsmap = IV_entries(dsmapIV) ;
/*
   ----------------------------------
   check compids[] and get the number 
   of domains and interface vertices
   ----------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
ndom = nPhi = 0 ;
for ( v = 0 ; v < nV ; v++ ) {
   if ( (d = compids[v]) < 0 ) {
      fprintf(stderr, 
              "\n fatal error in GPart_domSegMap(%p,%p,%p)"
              "\n compids[%d] = %d\n", gpart, pndom, pnseg,
              v, compids[v]) ;
      exit(-1) ;
   } else if ( d == 0 ) {
      nPhi++ ;
   } else if ( ndom < d ) {
      ndom = d ;
   }
}
*pndom = ndom ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n Inside GPart_domSegMap") ;
   fprintf(msgFile, "\n %d domains, %d Phi vertices", ndom, nPhi) ;
}
if ( ndom == 1 ) {
   IVfill(nV, dsmap, 0) ;
   *pndom = 1 ;
   *pnseg = 0 ;
   return(dsmapIV) ;
}
/*
   --------------------------------
   get the maps
   PhiToV : [0,nPhi) |---> [0,nV)
   VtoPhi : [0,nV)   |---> [0,nPhi)
   --------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
PhiToV = IVinit(nPhi, -1) ;
VtoPhi = IVinit(nV,   -1) ;
for ( v = 0, phi = 0 ; v < nV ; v++ ) {
   if ( (d = compids[v]) == 0 ) {
      PhiToV[phi] = v ;
      VtoPhi[v]   = phi++ ;
   }
}
if ( phi != nPhi ) {
   fprintf(stderr, 
           "\n fatal error in GPart_domSegMap(%p,%p,%p)"
           "\n phi = %d != %d = nPhi\n", 
           gpart, pndom, pnseg, phi, nPhi) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n PhiToV(%d) :", nPhi) ;
   IVfp80(msgFile, nPhi, PhiToV, 15, &ierr) ;
   fflush(msgFile) ;
}
if ( msglvl > 3 ) {
   fprintf(msgFile, "\n VtoPhi(%d) :", nV) ;
   IVfp80(msgFile, nV, VtoPhi, 15, &ierr) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------------
   create an IVL object, PhiByPowD, to hold lists from 
   the interface vertices to their adjacent domains
   ---------------------------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
dmark = IVinit(ndom+1, -1) ;
if ( nPhi >= ndom ) {
   list = IVinit(nPhi, -1) ;
} else {
   list = IVinit(ndom, -1) ;
}
PhiByPowD = IVL_new() ;
IVL_init1(PhiByPowD, IVL_CHUNKED, nPhi) ;
for ( phi = 0 ; phi < nPhi ; phi++ ) {
   v = PhiToV[phi] ;
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
/*
if ( phi == 0 ) {
   int   ierr ;
   fprintf(msgFile, "\n adj(%d,%d) = ", v, phi) ;
   IVfp80(msgFile, vsize, vadj, 15, &ierr) ;
   fflush(msgFile) ;
}
*/
   count = 0 ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      if ( (w = vadj[ii]) < nV  
        && (d = compids[w]) > 0 
        && dmark[d] != phi ) {
         dmark[d] = phi ;
         list[count++] = d ;
      }
   }
   if ( count > 0 ) {
      IVqsortUp(count, list) ;
      IVL_setList(PhiByPowD, phi, count, list) ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n PhiByPowD : interface x adjacent domains") ;
   IVL_writeForHumanEye(PhiByPowD, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------------------
   create an IVL object, PhiByPhi to hold lists
   from the interface vertices to interface vertices.
   (s,t) are in the list if (s,t) is an edge in the graph
   and s and t do not share an adjacent domain
   -------------------------------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
PhiByPhi = IVL_new() ;
IVL_init1(PhiByPhi, IVL_CHUNKED, nPhi) ;
offsets = IVinit(nPhi,  0)  ;
head    = IVinit(nPhi, -1) ;
link    = IVinit(nPhi, -1) ;
for ( phi1 = 0 ; phi1 < nPhi ; phi1++ ) {
   count = 0 ;
   if ( msglvl > 2 ) {
      v = PhiToV[phi1] ;
      Graph_adjAndSize(g, v, &vsize, &vadj) ;
      fprintf(msgFile, "\n checking out phi = %d, v = %d", phi1, v) ;
      fprintf(msgFile, "\n adj(%d) : ", v) ;
      IVfp80(msgFile, vsize, vadj, 10, &ierr) ;
   }
/*
   -------------------------------------------------------------
   get (phi1, phi0) edges that were previously put into the list
   -------------------------------------------------------------
*/
   if ( msglvl > 3 ) {
      if ( head[phi1] == -1 ) {
         fprintf(msgFile, "\n    no previous edges") ;
      } else {
         fprintf(msgFile, "\n    previous edges :") ;
      }
   }
   for ( phi0 = head[phi1] ; phi0 != -1 ; phi0 = nextphi ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, " %d", phi0) ;
      }
      nextphi = link[phi0] ;
      list[count++] = phi0 ;
      IVL_listAndSize(PhiByPhi, phi0, &size0, &adj0) ;
      if ( (ii = ++offsets[phi0]) < size0 ) {
/*
         ----------------------------
         link phi0 into the next list
         ----------------------------
*/
         phi2       = adj0[ii]   ;
         link[phi0] = head[phi2] ;
         head[phi2] = phi0       ;
      }
   }
/*
   --------------------------
   get new edges (phi1, phi2)
   --------------------------
*/
   IVL_listAndSize(PhiByPowD, phi1, &size1, &adj1) ;
   v = PhiToV[phi1] ;
   Graph_adjAndSize(g, v, &vsize, &vadj) ;
   for ( ii = 0 ; ii < vsize ; ii++ ) {
      if (  (w = vadj[ii]) < nV 
         && compids[w] == 0 
         && (phi2 = VtoPhi[w]) > phi1 ) {
         if ( msglvl > 3 ) {
            fprintf(msgFile, "\n    checking out phi2 = %d", phi2) ;
         }
/*
         --------------------------------------------------
         see if phi1 and phi2 have a common adjacent domain
         --------------------------------------------------
*/
         IVL_listAndSize(PhiByPowD, phi2, &size2, &adj2) ;
         adjdom = 0 ;
         jj1 = jj2 = 0 ;
         while ( jj1 < size1 && jj2 < size2 ) {
            if ( adj1[jj1] < adj2[jj2] ) {
               jj1++ ;
            } else if ( adj1[jj1] > adj2[jj2] ) {
               jj2++ ;
            } else {
               if ( msglvl > 3 ) {
                  fprintf(msgFile, ", common adj domain %d", adj1[jj1]) ;
               }
               adjdom = 1 ;
               break ;
            }
         }
         if ( adjdom == 0 ) {
            if ( msglvl > 3 ) {
               fprintf(msgFile, ", no adjacent domain") ;
            }
            list[count++] = phi2 ;
         }
      }
   }
   if ( count > 0 ) {
/*
      ---------------------
      set the list for phi1
      ---------------------
*/
      IVqsortUp(count, list) ;
      IVL_setList(PhiByPhi, phi1, count, list) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n    edge list for %d :", phi1) ;
         IVfp80(msgFile, count, list, 15, &ierr) ;
      }
      for ( ii = 0, phi2 = -1 ; ii < count ; ii++ ) {
         if ( list[ii] > phi1 ) {
            offsets[phi1] = ii ;
            phi2 = list[ii] ;
            break ;
         }
      }
      if ( phi2 != -1 ) {
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
                 "\n       linking %d into list for %d", phi1, phi2) ;
      }
         link[phi1] = head[phi2] ;
         head[phi2] = phi1       ;
      }
/*
      phi2 = list[0] ;
      link[phi1] = head[phi2] ;
      head[phi2] = phi1 ;
*/
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n PhiByPhi : interface x interface") ;
   IVL_writeForHumanEye(PhiByPhi, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------
   get the PhiToPsi map
   --------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
PhiToPsi = IVinit(nPhi, -1) ;
nPsi = 0 ;
for ( phi = 0 ; phi < nPhi ; phi++ ) {
   if ( PhiToPsi[phi] == -1 ) {
/*
      ---------------------------
      phi not yet mapped to a psi
      ---------------------------
*/
      first = last = 0 ;
      list[0] = phi ;
      PhiToPsi[phi] = nPsi ;
      while ( first <= last ) {
         phi2 = list[first++] ;
         IVL_listAndSize(PhiByPhi, phi2, &size, &adj) ;
         for ( ii = 0 ; ii < size ; ii++ ) {
            phi3 = adj[ii] ;
            if ( PhiToPsi[phi3] == -1 ) {
               PhiToPsi[phi3] = nPsi ;
               list[++last] = phi3 ; 
            }
         }
      }
      nPsi++ ;
   }
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n nPsi = %d", nPsi) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n PhiToPsi(%d) :", nPhi) ;
   IVfp80(msgFile, nPhi, PhiToPsi, 15, &ierr) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------
   create an IVL object, Psi --> 2^D
   ---------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
IVfill(nPsi, head, -1) ;
IVfill(nPhi, link, -1) ;
for ( phi = 0 ; phi < nPhi ; phi++ ) {
   psi = PhiToPsi[phi] ;
   link[phi] = head[psi] ;
   head[psi] =   phi     ;
}
PsiByPowD = IVL_new() ;
IVL_init1(PsiByPowD, IVL_CHUNKED, nPsi) ;
IVfill(ndom+1, dmark, -1) ;
for ( psi = 0 ; psi < nPsi ; psi++ ) {
   count = 0 ;
   for ( phi = head[psi] ; phi != -1 ; phi = link[phi] ) {
      v = PhiToV[phi] ;
      Graph_adjAndSize(g, v, &size, &adj) ;
      for ( ii = 0 ; ii < size ; ii++ ) {
         if (  (w = adj[ii]) < nV
            && (d = compids[w]) > 0 
            && dmark[d] != psi ) {
            dmark[d]      = psi ;
            list[count++] =  d  ;
         }
      }
   }
   if ( count > 0 ) {
      IVqsortUp(count, list) ;
      IVL_setList(PsiByPowD, psi, count, list) ;
   }
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n PsiByPowD(%d) :", nPhi) ;
   IVL_writeForHumanEye(PsiByPowD, msgFile) ;
   fflush(msgFile) ;
}
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
/*
   -------------------------------------
   now get the map Psi |---> Sigma that 
   is the equivalence map over PhiByPowD
   -------------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
PsiToSigma = IVL_equivMap1(PsiByPowD) ;
*pnseg = 1 + IVmax(nPsi, PsiToSigma, &ii) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n nSigma = %d", *pnseg) ;
   fprintf(msgFile, "\n PsiToSigma(%d) :", nPhi) ;
   IVfp80(msgFile, nPsi, PsiToSigma, 15, &ierr) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------------------
   now fill the map from the vertices to the domains and segments
   --------------------------------------------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
for ( v = 0 ; v < nV ; v++ ) {
   if ( (d = compids[v]) > 0 ) {
      dsmap[v] = d - 1 ;
   } else {
      phi      = VtoPhi[v] ;
      psi      = PhiToPsi[phi] ;
      sigma    = PsiToSigma[psi] ;
      dsmap[v] = ndom + sigma ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
IVL_free(PhiByPhi)  ;
IVL_free(PhiByPowD) ;
IVL_free(PsiByPowD) ;
IVfree(PhiToV)      ;
IVfree(VtoPhi)      ;
IVfree(dmark)       ;
IVfree(list)        ;
IVfree(PhiToPsi)    ;
IVfree(head)        ;
IVfree(link)        ;
IVfree(offsets)     ;
IVfree(PsiToSigma)  ;
icputimes++ ;
MARKTIME(cputimes[icputimes]) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n domain/segment map timings split") ;
   fprintf(msgFile, 
   "\n %9.5f : create the DSmap object"
   "\n %9.5f : get numbers of domain and interface vertices"
   "\n %9.5f : generate PhiToV and VtoPhi"
   "\n %9.5f : generate PhiByPowD"
   "\n %9.5f : generate PhiByPhi"
   "\n %9.5f : generate PhiToPsi"
   "\n %9.5f : generate PsiByPowD"
   "\n %9.5f : generate PsiToSigma"
   "\n %9.5f : generate dsmap"
   "\n %9.5f : free working storage"
   "\n %9.5f : total time",
   cputimes[1] - cputimes[0],
   cputimes[2] - cputimes[1],
   cputimes[3] - cputimes[2],
   cputimes[4] - cputimes[3],
   cputimes[5] - cputimes[4],
   cputimes[6] - cputimes[5],
   cputimes[7] - cputimes[6],
   cputimes[8] - cputimes[7],
   cputimes[9] - cputimes[8],
   cputimes[10] - cputimes[9],
   cputimes[11] - cputimes[0]) ;
}

return(dsmapIV) ; }

/*--------------------------------------------------------------------*/
