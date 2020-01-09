/*  testSemi.c  */

#include "../misc.h"
#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------
   get statistics for a semi-implicit solve

   created -- 97dec11, cca
   ----------------------------------------
*/
{
char     *inGraphFileName, *inETreeFileName, *inMapFileName ;
double   nA21, nL, nL11, nL22, nPhi, nV, t1, t2 ;
ETree    *etree ;
int      ii, inside, J, K, msglvl, nfront, nJ, 
         nvtx, rc, sizeJ, v, vsize, w ;
int      *adjJ, *frontmap, *map, *nodwghts, 
         *vadj, *vtxToFront, *vwghts ;
IV       *mapIV ;
IVL      *symbfacIVL ;
Graph    *graph ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 6 ) {
   fprintf(stdout, 
     "\n\n usage : %s msglvl msgFile GraphFile ETreeFile mapFile "
     "\n    msglvl    -- message level"
     "\n    msgFile   -- message file"
     "\n    GraphFile -- input graph file, must be *.graphf or *.graphb"
     "\n    ETreeFile -- input ETree file, must be *.etreef or *.etreeb"
     "\n    mapFile   -- input map IV file, must be *.ivf or *.ivb"
     "\n",
argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
inGraphFileName = argv[3] ;
inETreeFileName = argv[4] ;
inMapFileName   = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n GraphFile     -- %s" 
        "\n ETreeFile     -- %s" 
        "\n mapFile       -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inGraphFileName, inETreeFileName, inMapFileName) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
graph = Graph_new() ;
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
nvtx   = graph->nvtx ;
vwghts = graph->vwghts ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Graph object from file %s",
           inGraphFileName) ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
etree = ETree_new() ;
if ( strcmp(inETreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = ETree_readFromFile(etree, inETreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in etree from file %s",
        t2 - t1, inETreeFileName) ;
nfront     = ETree_nfront(etree) ;
tree       = ETree_tree(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
nL         = ETree_nFactorEntries(etree, 2) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, etree, inETreeFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading ETree object from file %s",
           inETreeFileName) ;
   ETree_writeForHumanEye(etree, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the map IV object
   -------------------------
*/
mapIV = IV_new() ;
if ( strcmp(inMapFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = IV_readFromFile(mapIV, inMapFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in mapIV from file %s",
        t2 - t1, inMapFileName) ;
map = IV_entries(mapIV) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, mapIV, inMapFileName) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading IV object from file %s",
           inMapFileName) ;
   IV_writeForHumanEye(mapIV, msgFile) ;
   fflush(msgFile) ;
}
nV = nPhi = 0 ;
if ( vwghts == NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      nV++ ;
      if ( map[v] == 0 ) {
         nPhi++ ;
      }
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      nV += vwghts[v] ;
      if ( map[v] == 0 ) {
         nPhi += vwghts[v] ;
      }
   }
}
fprintf(msgFile, "\n nPhi = %.0f, nV = %.0f", nPhi, nV) ;
/*
   -------------------------
   get the frontmap[] vector
   -------------------------
*/
frontmap = IVinit(nfront, -1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   J = vtxToFront[v] ;
   if ( frontmap[J] == -1 ) {
      frontmap[J] = map[v] ;
   } else if ( frontmap[J] != map[v] ) {
      fprintf(msgFile, "\n\n error, frontmap[%d] = %d, map[%d] = %d",
              J, frontmap[J], v, map[v]) ;
   }
}
/*
   ----------------------------------
   compute the symbolic factorization
   ----------------------------------
*/
symbfacIVL = SymbFac_initFromGraph(etree, graph) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n symbolic factorization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------
   compute the number of entries in L11 and L22
   --------------------------------------------
*/
nL11 = nL22 = 0 ;
for ( J = Tree_postOTfirst(tree) ;
      J != -1 ;
      J = Tree_postOTnext(tree, J) ) {
   nJ = nodwghts[J] ;
   if ( msglvl > 3 ) {
      fprintf(msgFile, "\n\n front %d, nJ = %d", J, nJ) ;
   }
   IVL_listAndSize(symbfacIVL, J, &sizeJ, &adjJ) ;
   for ( ii = 0, inside = 0 ; ii < sizeJ ; ii++ ) {
      w = adjJ[ii] ;
      K = vtxToFront[w] ;
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    w = %d, K = %d", w, K) ;
      }
      if ( K > J && frontmap[K] == frontmap[J] ) {
         inside += (vwghts == NULL) ? 1 : vwghts[w] ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, ", inside") ;
         }
      }
   }
   if ( frontmap[J] != 0 ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    inside = %d, adding %d to L11",
                 inside, nJ*nJ + 2*nJ*inside) ;
      }
      nL11 += nJ*nJ + 2*nJ*inside ;
   } else {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    inside = %d, adding %d to L22",
                 inside, nJ*nJ + 2*nJ*inside) ;
      }
      nL22 += nJ*nJ + 2*nJ*inside ;
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n |L| = %.0f, |L11| = %.0f, |L22| = %.0f", 
           nL, nL11, nL22) ;
}
/*
   ------------------------------------
   compute the number of entries in A21
   ------------------------------------
*/
nA21 = 0 ;
if ( vwghts != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      if ( map[v] == 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( map[v] != map[w] ) {
               if ( msglvl > 3 ) {
                  fprintf(msgFile, "\n A21 : v = %d, w = %d", v, w) ;
               }
               nA21 += vwghts[v] * vwghts[w] ;
            }
         }
      }
   }
} else {
   for ( v = 0 ; v < nvtx ; v++ ) {
      if ( map[v] == 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( map[v] != map[w] ) {
               if ( msglvl > 3 ) {
                  fprintf(msgFile, "\n A21 : v = %d, w = %d", v, w) ;
               }
               nA21++ ;
            }
         }
      }
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n |L| = %.0f, |L11| = %.0f, |L22| = %.0f, |A21| = %.0f", 
           nL, nL11, nL22, nA21) ;
   fprintf(msgFile, 
      "\n storage: explicit = %.0f, semi-implicit = %.0f, ratio = %.3f"
      "\n opcount: explicit = %.0f, semi-implicit = %.0f, ratio = %.3f",
      nL, nL11 + nA21 + nL22, 
      nL/(nL11 + nA21 + nL22),
      2*nL, 4*nL11 + 2*nA21 + 2*nL22,
      2*nL/(4*nL11 + 2*nA21 + 2*nL22)) ;
   fprintf(msgFile, "\n ratios %8.3f %8.3f %8.3f",
           nPhi/nV,
           nL/(nL11 + nA21 + nL22),
           2*nL/(4*nL11 + 2*nA21 + 2*nL22)) ;
}
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
ETree_free(etree) ;
IV_free(mapIV) ;
IVL_free(symbfacIVL) ;
IVfree(frontmap) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
