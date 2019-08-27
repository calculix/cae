/*  RBviaDDsep.c  */

#include "../DDsepInfo.h"
#include "../GPart.h"
#include "../../timings.h"

#define NSPLIT 8 

#define BE_CAUTIOUS 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   prototype for static visit() method
   -----------------------------------
*/
static void visit ( GPart *gpart, int map[], int vpar[], 
                    IV *DDmapIV, DDsepInfo *info) ;

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   this method constructs recursive partition of a graph.
   it returns a DSTree object to represent the splitting of
   the tree into subgraphs.
   the info object contains the information needed by the
   DDsep algorithm.

   created -- 96feb24, cca
   --------------------------------------------------------
*/
DSTree *
GPart_RBviaDDsep (
   GPart       *gpart,
   DDsepInfo   *info
) {
double   t0, t1, t2, t3 ;
DSTree   *dstree ;
FILE     *msgFile ;
GPart    *child ;
int      ierr, msglvl, nvtx ;
int      *map, *vpar ;
IV       *DDmapIV, *DSmapIV ;
Tree     *tree ;
/*
   ---------------
   check the input
   ---------------
*/
MARKTIME(t0) ;
if (    gpart == NULL 
     || (nvtx = gpart->nvtx) <= 0 
     || info == NULL
   ) {
   fprintf(stderr, "\n fatal error in GPart_RBviaDDsep(%p,%p)"
           "\n bad input\n", gpart, info) ;
   exit(-1) ;
}
msglvl  = gpart->msglvl  ;
msgFile = gpart->msgFile ;
/*
   -------------------------
   check that gpart is a root
   -------------------------
*/
if ( gpart->par != NULL ) {
   fprintf(stderr, "\n fatal error in GPart_RBviaDDsep(%p,%p)"
           "\n gpart must be a root \n", gpart, info) ;
   exit(-1) ;
}
/*
   -----------------------------
   create map and parent vectors
   -----------------------------
*/
vpar = IVinit(nvtx, -1) ;
DSmapIV = IV_new() ;
IV_init(DSmapIV, nvtx, NULL) ;
map = IV_entries(DSmapIV) ;
IVfill(nvtx, map, -1) ;
info->ntreeobj = 0 ;

if ( info->DDoption == 2 ) {
/*
   ------------------------------------------------
   get one global domain decomposition via fishnet,
   then project it on each subgraph
   ------------------------------------------------
*/
   MARKTIME(t1) ;
   GPart_DDviaFishnet(gpart, info->freeze, info->minweight,
                      info->maxweight, info->seed) ;
   DDmapIV = IV_new() ;
   IV_init(DDmapIV, nvtx, NULL) ;
   IV_copy(DDmapIV, &gpart->compidsIV) ;
   IV_fill(&gpart->compidsIV, 1) ;
   MARKTIME(t2) ;
   info->cpuDD += t2 - t1 ;
} else {
   DDmapIV = NULL ;
}
/*
   ---------------------------------------------
   split the graph into its connected components
   ---------------------------------------------
*/
MARKTIME(t1) ;
GPart_split(gpart) ;
MARKTIME(t2) ;
info->cpuSplit += t2 - t1 ;
if ( msglvl > 2 && msgFile != NULL ) {
   fprintf(msgFile, "\n after initial split, ncomp = %d",
           gpart->ncomp) ;
   fflush(msgFile) ;
}
if ( gpart->ncomp > 0 ) {
   for ( child = gpart->fch ; child != NULL ; child = child->sib ) {
      child->id = info->ntreeobj++ ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n\n ### component %d", child->id) ;
         Graph_writeStats(child->g, msgFile) ;
         if ( msglvl > 3 && msgFile != NULL ) {
            Graph_writeForHumanEye(child->g, msgFile) ;
            if ( IV_size(&child->vtxMapIV) > 0 ) {
               fprintf(msgFile, 
                       "\n vtxMap(%d) :", child->nvtx) ;
               IV_fp80(&child->vtxMapIV, msgFile, 80, &ierr) ;
            }
         }
      }
      fflush(msgFile) ;
   }
}
/*
   ------------------------------
   visit each connected component
   ------------------------------
*/
if ( gpart->fch != NULL ) {
   while ( (child = gpart->fch) != NULL ) {
      gpart->fch = child->sib ;
      visit(child, map, vpar, DDmapIV, info) ;
      Graph_free(child->g) ;
      GPart_free(child) ;
   }
} else {
   gpart->id = info->ntreeobj++ ;
   visit(gpart, map, vpar, DDmapIV, info) ;
}
/*
   ------------------------
   create the DSTree object
   ------------------------
*/
tree = Tree_new() ;
Tree_init2(tree, info->ntreeobj, vpar) ;
dstree = DSTree_new() ;
DSTree_init2(dstree, tree, DSmapIV) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
IVfree(vpar) ;
MARKTIME(t3) ;
info->cpuTotal = t3 - t0 ;

return(dstree) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   this method is called for each connected subgraph in the
   subgraph tree.

   map  -- vector mapping vertices to subgraph number
   vpar -- vector mapping a subgraph to its parent subgraph
   info -- information structure for DDsep

   created -- 96feb24, cca
   --------------------------------------------------------
*/
static void
visit ( 
   GPart       *gpart,
   int         map[], 
   int         vpar[],
   IV          *DDmapIV,
   DDsepInfo   *info
) {
double   t1, t2 ;
double   cpus[3] ;
FILE     *msgFile ;
GPart    *child, *par ;
Graph    *g ;
int      ierr, ii, maxweight, minweight, msglvl,
         nvbnd, nvtot, nvtx, parnvtot, totvwght, v, vsize, width ;
int      *compids, *parmap, *vadj, *vtxMap ;
/*
   -------------------------------
   extract dimensions and pointers
   -------------------------------
*/
g       = gpart->g ;
nvtx    = g->nvtx ;
nvbnd   = g->nvbnd ;
nvtot   = nvtx + nvbnd  ;
compids = IV_entries(&gpart->compidsIV) ;
vtxMap  = IV_entries(&gpart->vtxMapIV)  ;
msgFile = gpart->msgFile ;
msglvl  = gpart->msglvl  ;
/*
   ------------------------------
   compute the total graph weight
   ------------------------------
*/
if ( g->type % 2 == 0 ) {
   totvwght = nvtx ;
} else {
   totvwght = IVsum(nvtx, g->vwghts) ;
}
/*
   -----------------
   optional messages
   -----------------
*/
if ( msglvl > 1 && msgFile != NULL ) {
   fprintf(msgFile, 
           "\n\n ### inside visit(%d), parent = %d"
             "\n     nvtx = %d, nvbnd = %d, nvtot = %d, totvwght = %d",
           gpart->id, (gpart->par != NULL) ? gpart->par->id : -1,
           nvtx, nvbnd, nvtot, totvwght) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 && msgFile != NULL ) {
   Graph_writeForHumanEye(gpart->g, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   set the vertex map to be global
   -------------------------------
*/
if (  (par = gpart->par) != NULL 
   && (parmap = IV_entries(&par->vtxMapIV)) != NULL ) {
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n before changing map") ;
      fprintf(msgFile, "\n vtxMapIV") ;
      IV_writeForHumanEye(&gpart->vtxMapIV, msgFile) ;
      fprintf(msgFile, "\n parVtxMapIV") ;
      IV_writeForHumanEye(&par->vtxMapIV, msgFile) ;
      fflush(msgFile) ;
      IV_fp80(&gpart->vtxMapIV, msgFile, 80, &ierr) ;
      fflush(msgFile) ;
   }
   parnvtot = par->nvtx + par->nvbnd ;
   for ( v = 0 ; v < nvtot ; v++ ) {
#if BE_CAUTIOUS > 0
      if ( vtxMap[v] < 0 || vtxMap[v] >= parnvtot ) {
         fprintf(stderr, 
              "\n error changing map, vtxMap[[%d] = %d, parnvtot = %d",
              v, vtxMap[v], parnvtot) ;
         exit(-1) ;
      }
#endif
      vtxMap[v] = parmap[vtxMap[v]] ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n after changing map") ;
      IV_fp80(&gpart->vtxMapIV, msgFile, 80, &ierr) ;
      fflush(msgFile) ;
   }
}
if ( totvwght <= info->maxcompweight ) {
   gpart->ncomp = 1 ;
} else {
/*
   ----------------------------------------------
   try to find a bisector via the DDsep algorithm
   ----------------------------------------------
*/
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n try to find a bisector") ;
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   switch ( info->DDoption ) {
   case 1 :
/*
      ----------------------------------------------------------
      use the fishnet algorithm to find the domain decomposition
      ----------------------------------------------------------
*/
      if ( info->maxweight * NSPLIT <= totvwght ) {
         minweight = info->minweight ;
         maxweight = info->maxweight ;
      } else {
         maxweight = totvwght / NSPLIT ;
         if ( maxweight < 2 ) {
            maxweight = 2 ;
         }
         minweight = maxweight / 2 ;
      }
      if ( msglvl > 2 ) {
         fprintf(msgFile, 
        "\n calling DDviaFishnet with minweight = %d, maxweight = %d",
                 minweight, maxweight) ;
         fflush(msgFile) ;
      }
      GPart_DDviaFishnet(gpart, info->freeze, minweight, 
                         maxweight, info->seed) ;
      if ( msglvl > 2 ) {
         fprintf(msgFile, "\n return from DDviaFishnet") ;
         fflush(msgFile) ;
      }
      break ;
   case 2 :
/*
      -----------------------------------------------------------------
      use projection from a global map to find the domain decomposition
      -----------------------------------------------------------------
*/
      GPart_DDviaProjection(gpart, DDmapIV) ;
      break ;
   }
   MARKTIME(t2) ;
   info->cpuDD += t2 - t1 ;
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n after DD: %d domains", gpart->ncomp) ;
   }
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile, "\n partition weights :") ;
      IV_fp80(&gpart->cweightsIV, msgFile, 25, &ierr) ;
   }
#if BE_CAUTIOUS > 0
   if ( GPart_validVtxSep(gpart) == 1 ) {
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(stdout, "\n DD: valid vertex separator ") ;
      }
   } else {
      fprintf(stderr, "\n DD: invalid vertex separator ") ;
      exit(-1) ;
   }
#endif
   if ( gpart->ncomp > 1 ) {
/*
      --------------------------------
      find an initial bisector via BKL
      --------------------------------
*/
      GPart_TwoSetViaBKL(gpart, info->alpha, info->seed, cpus) ;
      info->cpuMap += cpus[0] ;
      info->cpuBPG += cpus[1] ;
      info->cpuBKL += cpus[2] ;
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n BKL final weights   : ") ;
         IV_fp80(&gpart->cweightsIV, msgFile, 25, &ierr) ;
      }
#if BE_CAUTIOUS > 0
      if ( GPart_validVtxSep(gpart) == 1 ) {
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(stdout, "\n BKL: valid vertex separator ") ;
         }
      } else {
         fprintf(stderr, "\n BKL: invalid vertex separator ") ;
         exit(-1) ;
      }
#endif
   }
   if ( gpart->ncomp > 1 ) {
/*
      ---------------------------------------
      find an improved bisector via smoothing
      ---------------------------------------
*/
      MARKTIME(t1) ;
      if ( info->nlayer <= 2 ) {
         GPart_smoothBy2layers(gpart, info->nlayer, info->alpha) ;
      } else if ( (width = info->nlayer/2) > 0 ) {
         GPart_smoothBisector(gpart, width, info->alpha) ;
      }
      MARKTIME(t2) ;
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n smoothed weights          : ") ;
         IV_fp80(&gpart->cweightsIV, msgFile, 25, &ierr) ;
      }
#if BE_CAUTIOUS > 0
      if ( GPart_validVtxSep(gpart) == 1 ) {
         if ( msglvl > 2 && msgFile != NULL ) {
            fprintf(stdout, "\n smoothed: valid vertex separator ") ;
         }
      } else {
         fprintf(stderr, "\n smoothed: invalid vertex separator ") ;
         exit(-1) ;
      }
#endif
      info->cpuSmooth += t2 - t1 ;
   } 
   if ( gpart->ncomp > 1 ) {
/*
      ----------------------------------
      split the subgraph into components
      ----------------------------------
*/
      MARKTIME(t1) ;
      GPart_split(gpart) ;
      MARKTIME(t2) ;
      info->cpuSplit += t2 - t1 ;
      if ( msglvl > 1 && msgFile != NULL ) {
         fprintf(msgFile, "\n SPLIT weights       : ") ;
         IV_fp80(&gpart->cweightsIV, msgFile, 20, &ierr) ;
         fflush(msgFile) ;
      }
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile, "\n compids") ;
         IV_fp80(&gpart->compidsIV, msgFile, 80, &ierr) ;
         fflush(msgFile) ;
      }
   }
}
if ( gpart->ncomp > 1 ) {
/*
   ---------------------------------------
   a separator was found and the graph has 
   been split into connected components
   ---------------------------------------
*/
   for ( child = gpart->fch ; child != NULL ; child = child->sib) {
      child->id = info->ntreeobj++ ;
      vpar[child->id] = gpart->id ;
   }
   if ( msglvl > 2 && msgFile != NULL ) {
      fprintf(msgFile,
              "\n after initial split, ncomp = %d", gpart->ncomp) ;
      for ( child = gpart->fch ; child != NULL ; child = child->sib) {
         fprintf(msgFile, "\n\n ### component %d", child->id) ;
         Graph_writeStats(child->g, msgFile) ;
         if ( msglvl > 3 && msgFile != NULL ) {
            Graph_writeForHumanEye(child->g, msgFile) ;
            if ( IV_size(&child->vtxMapIV) > 0 ) {
               fprintf(msgFile, 
                       "\n vtxMap(%d) :", child->nvtx) ;
               IV_fp80(&child->vtxMapIV, msgFile,  80, &ierr) ;
            }
         }
      }
      fflush(msgFile) ;
   }
/*
   ------------------------------
   visit each connected component
   ------------------------------
*/
   while ( (child = gpart->fch) != NULL ) {
      gpart->fch = child->sib ;
      visit(child, map, vpar, DDmapIV, info) ;
      if ( msglvl > 2 && msgFile != NULL ) {
         fprintf(msgFile,
                 "\n return from visiting child %d", child->id) ;
         fflush(msgFile) ;
      }
      Graph_free(child->g) ;
      GPart_free(child) ;
   }
/*
   -----------------------------------------------------------
   map adjacency structure of the separator back to the global
   numbering and set the map for the separator vertices
   -----------------------------------------------------------
*/
   if ( gpart->par != NULL ) {
      for ( v = 0 ; v < nvtx ; v++ ) {
         if ( compids[v] == 0 ) {
            Graph_adjAndSize(g, v, &vsize, &vadj) ;
            for ( ii = 0 ; ii < vsize ; ii++ ) {
#if BE_CAUTIOUS > 0
            if ( vadj[ii] < 0 || vadj[ii] >= nvtot ) {
               fprintf(stderr, 
                       "\n 1.0 whoa, error, vadj[[%d] = %d, nvtot = %d",
                       ii, vadj[ii], nvtot) ;
               exit(-1) ;
            }
#endif
               vadj[ii] = vtxMap[vadj[ii]] ;
            }
            map[vtxMap[v]] = gpart->id ;
         }
      }
   } else {
      for ( v = 0 ; v < nvtx ; v++ ) {
         if ( compids[v] == 0 ) {
            map[v] = gpart->id ;
         }
      }
   }
} else {
/*
   -----------------------------------------------------
   the subgraph was not split. map the adjacency back to
   the global numbering and set the map for the vertices
   -----------------------------------------------------
*/
   if ( msglvl > 1 && msgFile != NULL ) {
      fprintf(msgFile, "\n this subgraph is a domain") ;
   }
   if ( gpart->par != NULL ) {
      for ( v = 0 ; v < nvtx ; v++ ) {
         Graph_adjAndSize(g, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
#if BE_CAUTIOUS > 0
            if ( vadj[ii] < 0 || vadj[ii] >= nvtot ) {
               fprintf(stderr, 
                       "\n 2.0 whoa, error, vadj[[%d] = %d, nvtot = %d",
                       ii, vadj[ii], nvtot) ;
               exit(-1) ;
            }
#endif
            vadj[ii] = vtxMap[vadj[ii]] ;
         }
      }
      for ( v = 0 ; v < nvtx ; v++ ) {
         map[vtxMap[v]] = gpart->id ;
      }
   } else {
      for ( v = 0 ; v < nvtx ; v++ ) {
         map[v] = gpart->id ;
      }
   }
}

return ; }

/*--------------------------------------------------------------------*/
