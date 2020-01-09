/*  DDviaFishnet.c  */

#include "../GPart.h"
#include "../../IIheap.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   declarations for static functions
   ---------------------------------
*/
static int
GPart_freeze (
   GPart    *gpart,
   double   freeze,
   int      extdegs[]
) ;
static void
GPart_indpSetGrowth (
   GPart   *gpart,
   int     maxWeight,
   int     seed
) ;
static void
GPart_absDomains (
   GPart   *gpart,
   int     minweight
) ;
static void
GPart_absBoundary (
   GPart   *gpart
) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to construct and return a multisector
              using the fishnet algorithm

   freeze    -- parameter used to keep nodes of high degree
                in the interface
   minweight -- minimum weight for a domain
   maxweight -- maximum weight for a domain
   seed      -- random number seed

   created -- 96feb16, cca
   ---------------------------------------------------------------
*/
void
GPart_DDviaFishnet (
   GPart    *gpart,
   double   freeze,
   int      minweight,
   int      maxweight,
   int      seed
) {
double   cpus[6], t1, t2 ;
int      nV, v ;
int      *extdegs ;
/*
   ---------------
   check the input
   ---------------
*/
if (  gpart == NULL || freeze < 0.0 
   || minweight < 0 || maxweight < 0 || minweight >= maxweight ) {
   fprintf(stderr, 
          "\n fatal error in GPart_DDviaFishnet(%p,%f,%d,%d,%d)"
          "\n bad input\n", gpart, freeze, minweight, maxweight, seed) ;
   exit(-1) ;
}
/*
   ----------------------------
   compute the external degrees
   ----------------------------
*/
MARKTIME(t1) ;
nV = gpart->g->nvtx ;
extdegs = IVinit(nV, 0) ;
for ( v = 0 ; v < nV ; v++ ) {
   extdegs[v] = Graph_externalDegree(gpart->g, v) ;
}
MARKTIME(t2) ;
cpus[0] = t2 - t1 ;
/*
   -------------------
   freeze the vertices
   -------------------
*/
MARKTIME(t1) ;
GPart_freeze(gpart, freeze, extdegs) ;
MARKTIME(t2) ;
cpus[1] = t2 - t1 ;
/*
   ---------------------------------------------
   grow the partition via independent set growth
   ---------------------------------------------
*/
MARKTIME(t1) ;
GPart_indpSetGrowth(gpart, maxweight, seed) ;
IVfree(extdegs) ;
MARKTIME(t2) ;
cpus[2] = t2 - t1 ;
if ( gpart->ncomp == 1 ) {
   IV_fill(&gpart->compidsIV, 1) ;
   return ;
}
/*
   ------------------------
   absorb the small domains
   ------------------------
*/
MARKTIME(t1) ;
GPart_absDomains(gpart, minweight) ;
MARKTIME(t2) ;
cpus[3] = t2 - t1 ;
if ( gpart->ncomp <= 1 ) {
   IV_fill(&gpart->compidsIV, 1) ;
   return  ;
}
/*
   --------------------------
   absorb the excess boundary
   --------------------------
*/
MARKTIME(t1) ;
GPart_absBoundary(gpart) ;
MARKTIME(t2) ;
cpus[4] = t2 - t1 ;

if ( gpart->msglvl > 1 ) {
   fprintf(gpart->msgFile, 
           "\n FISHNET CPU breakdown"
           "\n CPU %8.3f : compute external degrees"
           "\n CPU %8.3f : freeze vertices of high degree"
           "\n CPU %8.3f : independent set growth"
           "\n CPU %8.3f : absorb small domains"
           "\n CPU %8.3f : absorb excess boundary",
           cpus[0], cpus[1], cpus[2], cpus[3], cpus[4]) ;
}

return ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------
   freeze vertices into the interface that 
   have external degree >= freeze * cutoff

   return value -- # of frozen vertices
 
   created -- 95oct05, cca
   ---------------------------------------
*/
static int
GPart_freeze (
   GPart    *gpart,
   double   freeze,
   int      extdegs[]
) {
Graph   *g ;
int     cutoff, iv, median, nfrozen, nvtx ;
int     *compids, *vids ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL || extdegs == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_freeze(%p,%f,%p)"
           "\n bad input\n", gpart, freeze, extdegs) ;
   exit(-1) ;
 }
nvtx = gpart->nvtx ;
compids = IV_entries(&gpart->compidsIV) ;
/*
   -------------------------------------------------------
   sort the vertices in ascending order of external degree
   -------------------------------------------------------
*/
vids = IVinit(nvtx, 0) ;
IVramp(nvtx, vids, 0, 1) ;
if ( gpart->msglvl > 3 ) {
   int v ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      fprintf(gpart->msgFile, 
              "\n vertex %d, external degree %d", v, extdegs[v]) ;
      fflush(gpart->msgFile) ;
   }
}
IV2qsortUp(nvtx, extdegs, vids) ;
/*
   --------------------------------
   get the median and cutoff values
   --------------------------------
*/
median = extdegs[nvtx/2] ;
cutoff = (int) (freeze * median) ;
if ( gpart->msglvl > 2 ) {
   fprintf(gpart->msgFile, "\n median = %d, cutoff = %d", median, cutoff) ;
   fflush(gpart->msgFile) ;
}
/*
   -------------------------------------------------------
   freeze vertices whose degree is greater than the cutoff
   -------------------------------------------------------
*/
nfrozen = 0 ;
for ( iv = nvtx - 1 ; iv >= 0 ; iv-- ) {
   if ( extdegs[iv] < cutoff ) {
      break ;
   }
   compids[vids[iv]] = 0 ;
   nfrozen++ ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(vids)    ;

return(nfrozen) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------------
   set up the partition via independent set growth

   maxWeight -- maximum weight for any domain
   seed      -- random number seed, if seed > 0 then
                the vertices are visited in some random order,
                else the vertices are visited 0, 1, ..., nvtx - 1

   note : the gpart->compids[] vector is assumed to be initialized
          outside this method. any vertex v such that compids[v] == -1
          is eligible to be put into a domain. this allows the user
          to pre-specify vertices to be in the separator (by setting
          compids[v] = -1 prior to calling this method), e.g.,
          vertices of abnormally high degree should be forced to be
          separator vertices.
  

   created -- 95oct05, cca
   -------------------------------------------------------------------
*/
static void
GPart_indpSetGrowth (
   GPart   *gpart,
   int     maxWeight,
   int     seed
) {
Graph   *g ;
int     domweight, i, iv, last, ndom, now, nvtx, v, vsize, w ;
int     *compids, *cweights, *list, *vadj, *vids, *vwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL || maxWeight < 0 ) {
   fprintf(stderr, "\n fatal error in GPart_indpSepGrowth(%p,%d,%d)"
           "\n bad input\n", gpart, maxWeight, seed) ;
   exit(-1) ;
}
vwghts = g->vwghts ;
/*
   --------------------------------
   set all component ids != 0 to -1
   --------------------------------
*/
nvtx    = gpart->nvtx ;
compids = IV_entries(&gpart->compidsIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] != 0 ) {
      compids[v] = -1 ;
   }
}
/*
   -----------------------------------------------------------------
   set up the vector that determines the order to visit the vertices
   -----------------------------------------------------------------
*/
vids = IVinit2(nvtx) ;
IVramp(nvtx, vids, 0, 1) ;
if ( seed > 0 ) {
   IVshuffle(nvtx, vids, seed) ;
}
/*
   ------------------------------
   set up the working list vector
   ------------------------------
*/
list = IVinit(nvtx, -1) ;
/*
   -----------------------------------
   visit the vertices and grow domains
   -----------------------------------
*/
ndom = 0 ;
for ( iv = 0 ; iv < nvtx ; iv++ ) {
   v = vids[iv] ;
   if ( gpart->msglvl > 4 ) {
      fprintf(gpart->msgFile, "\n\n visiting v = %d, compids[%d] = %d", 
              v, v, compids[v]) ;
   }
   if ( compids[v] == -1 ) {
/*
      -----------------------------------------
      eligible vertex has not yet been visited,
      create a new domain with v as the seed
      -----------------------------------------
*/
      ndom++ ;
      if ( gpart->msglvl > 3 ) {
         fprintf(gpart->msgFile, 
                 "\n\n domain %d : seed vertex %d", ndom, v) ;
         fflush(gpart->msgFile) ;
      }
      domweight = 0 ;
      now = last = 0 ;
      list[0] = v ;
      while ( now <= last ) {
         v = list[now++] ;
         if ( gpart->msglvl > 4 ) {
            fprintf(gpart->msgFile, 
                    "\n    adding %d to domain %d, weight %d",
                    v, ndom, ((vwghts != NULL) ? vwghts[v] : 1)) ;
            fflush(gpart->msgFile) ;
         }
/*
         ---------------------
         add v to domain icomp
         ---------------------
*/
         compids[v] = ndom ;
         domweight += (vwghts != NULL) ? vwghts[v] : 1 ;
         Graph_adjAndSize(g, v, &vsize, &vadj) ;
         for ( i = 0 ; i < vsize ; i++ ) {
            if ( (w = vadj[i]) < nvtx && compids[w] == -1 ) {
/*
               -----------------------------------------
               w has not yet been touched, put on list
               set compids[w] = -2 to make sure it isn't
               put on the list twice
               -----------------------------------------
*/
               compids[w]   = -2 ;
               list[++last] =  w ;
            }
         }
         if ( domweight >= maxWeight ) {
/*
            ---------------------------------------------
            domain is large enough, mark all the rest of 
            the vertices in the list as boundary vertices
            ---------------------------------------------
*/
            while ( now <= last ) {
               w = list[now++] ;
               if ( gpart->msglvl > 4 ) {
                  fprintf(gpart->msgFile, 
                          "\n    adding %d to interface, weight %d",
                          w, ((vwghts != NULL) ? vwghts[w] : 1)) ;
                  fflush(gpart->msgFile) ;
               }
               compids[w] = 0 ;
            }
         }
      }
      if ( gpart->msglvl > 2 ) {
         fprintf(gpart->msgFile, 
                 "\n domain %d, weight %d", ndom, domweight) ;
         fflush(gpart->msgFile) ;
      }
   }
}
/*
   ---------------------------------------------------------------------
   set the number of components and resize the component weights vector
   ---------------------------------------------------------------------
*/
gpart->ncomp = ndom ;
IV_setSize(&gpart->cweightsIV, 1 + ndom) ;
IV_fill(&gpart->cweightsIV, 0) ;
cweights = IV_entries(&gpart->cweightsIV) ;
/*
   -----------------------------------------------
   set the weights of the interface and components
   -----------------------------------------------
*/
if ( vwghts != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      cweights[compids[v]] += vwghts[v] ;
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      cweights[compids[v]]++ ;
   }
}
/*
   -------------------------
   free the two work vectors
   -------------------------
*/
IVfree(list) ;
IVfree(vids) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   absorb small domains into the interface.
   this method assumes that the component weights vector has been set

   95oct05, cca
   ------------------------------------------------------------------
*/
static void
GPart_absDomains (
   GPart   *gpart,
   int     minweight
) {
Graph   *g ;
int     c, ierr, ndom, nnewdom, nvtx, v ;
int     *compids, *cweights, *dmap, *head, *link, *vwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_absDomains(%p,%d)"
           "\n bad input\n", gpart, minweight) ;
   exit(-1) ;
}
nvtx     = gpart->nvtx  ;
vwghts   = g->vwghts    ;
ndom     = gpart->ncomp ;
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
/*
   ----------------------------------------------
   get the linked list for vertices and component
   ----------------------------------------------
*/
head = IVinit(ndom+1, -1) ;
link = IVinit(nvtx,   -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   c       = compids[v] ;
   link[v] = head[c]    ;
   head[c] = v          ;
}
dmap    = IVinit(ndom+1, -1) ;
nnewdom = 0 ;
dmap[0] = 0 ;
for ( c = 1 ; c <= ndom ; c++ ) {
   if ( cweights[c] < minweight ) {
      if ( gpart->msglvl > 2 ) {
         fprintf(gpart->msgFile, 
                 "\n interface absorbs component %d, weight %d",
                 c, cweights[c]) ;
         fflush(gpart->msgFile) ;
      }
      for ( v = head[c] ; v != -1 ; v = link[v] ) {
         compids[v] = 0 ;
      }
      cweights[0] += cweights[c] ;
      cweights[c] =  0 ;
      dmap[c]     =  0 ;
   } else {
      dmap[c] = ++nnewdom ;
   }
   if ( gpart->msglvl > 2 ) {
      fprintf(gpart->msgFile, "\n dmap[%d] = %d", c, dmap[c]) ;
      fflush(gpart->msgFile) ;
   }
}
if ( nnewdom != ndom ) {
/*
   ------------------------
   set the new # of domains
   ------------------------
*/
   gpart->ncomp = nnewdom ;
/*
   -------------------------
   set the new component ids
   -------------------------
*/
   if ( gpart->msglvl > 3 ) {
      fprintf(gpart->msgFile, "\n old component ids") ;
      IVfp80(gpart->msgFile, nvtx, compids, 80, &ierr) ;
      fflush(gpart->msgFile) ;
   }
   for ( v = 0 ; v < nvtx ; v++ ) {
      c = compids[v] ;
      compids[v] = dmap[c] ;
   }
   if ( gpart->msglvl > 3 ) {
      fprintf(gpart->msgFile, "\n new component ids") ;
      IVfp80(gpart->msgFile, nvtx, compids, 80, &ierr) ;
      fflush(gpart->msgFile) ;
   }
/*
   -----------------------------
   set the new component weights
   -----------------------------
*/
   if ( gpart->msglvl > 2 ) {
      fprintf(gpart->msgFile, "\n old cweights") ;
      IVfp80(gpart->msgFile, ndom+1, cweights, 80, &ierr) ;
      fflush(gpart->msgFile) ;
   }
   for ( c = 1 ; c <= ndom ; c++ ) {
      if ( dmap[c] != 0 ) {
         cweights[dmap[c]] = cweights[c] ;
      }
   }
   IV_setSize(&gpart->cweightsIV, nnewdom) ;
   if ( gpart->msglvl > 2 ) {
      fprintf(gpart->msgFile, "\n new cweights") ;
      IVfp80(gpart->msgFile, nnewdom+1, cweights, 80, &ierr) ;
      fflush(gpart->msgFile) ;
   }
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(head) ;
IVfree(link) ;
IVfree(dmap) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   absorb the excess boundary

   created -- 95oct05, cca
   revised -- 95oct19, cca
      very simple scheme
   ------------------------------------------------------------
*/
static void
GPart_absBoundary (
   GPart   *gpart
) {
Graph    *g ;
int      count, domid, ierr, ii, oldcount, ntest, nvtx, rc, v, vwght ;
int      *compids, *cweights, *list, *vwghts ;
/*
   ---------------
   check the input
   ---------------
*/
if ( gpart == NULL || (g = gpart->g) == NULL ) {
   fprintf(stderr, "\n fatal error in GPart_absBoundary(%p)"
           "\n bad input\n", gpart) ;
   exit(-1) ;
}
nvtx     = gpart->nvtx ;
compids  = IV_entries(&gpart->compidsIV)  ;
cweights = IV_entries(&gpart->cweightsIV) ;
vwghts   = gpart->g->vwghts ;
/*
   ----------------------
   create working storage
   ----------------------
*/
list = IVinit(nvtx, -1) ;
/*
   -----------------------------------------
   load all interface vertices into the list
   ------------------------------------------
*/
count = 0 ;
for ( v = 0 ; v < nvtx ; v++ ) {
   if ( compids[v] == 0 ) {
      list[count++] = v ;
   }
}
oldcount = -1 ;
while ( count > 0 ) {
   if ( gpart->msglvl > 2 ) {
      fprintf(gpart->msgFile, "\n\n new pass, count = %d", count) ;
   }
   ntest = count ;
   count = 0 ;
   for ( ii = 0 ; ii < ntest ; ii++ ) { 
      v = list[ii] ;
      rc = GPart_vtxIsAdjToOneDomain(gpart, v, &domid) ;
      if ( rc == 1 ) {
/*
         -----------------------------------------
         transfer interface vertex v into domain c
         -----------------------------------------
*/
         compids[v] = domid ;
         vwght = (vwghts != NULL) ? vwghts[v] : 1 ;
         cweights[0]     -= vwght ;
         cweights[domid] += vwght ;
         if ( gpart->msglvl > 3 ) {
            fprintf(gpart->msgFile, 
                 "\n    moving vertex %d with weight %d to domain %d"
                 "\n    now, cweights[0] = %d, cweights[%d] = %d",
                 v, vwght, domid, cweights[0], domid, cweights[domid]) ;
            fflush(gpart->msgFile) ;
         }
      } else if ( domid == -1 ) {
/*
         ----------------------------------------------------------
         vertex v is still not adjacent to any domain, keep on list
         ----------------------------------------------------------
*/
         if ( gpart->msglvl > 3 ) {
            fprintf(gpart->msgFile, 
                    "\n    keeping vertex %d on list", v) ;
         }
         list[count++] = v ;
      }
   }
   if ( count == oldcount ) {
      break ;
   } else {
      oldcount = count ;
   }
}
IVfree(list) ; 

return ; }
         
/*--------------------------------------------------------------------*/
