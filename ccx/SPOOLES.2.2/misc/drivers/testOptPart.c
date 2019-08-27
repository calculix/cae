/*  testOptPart.c  */

#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------
   (1) read in an ETree object.
   (2) read in an Graph object.
   (3) find the optimal domain/schur complement partition
       for a semi-implicit factorization
   
   created -- 96oct03, cca
   ------------------------------------------------------
*/
{
char     *inETreeFileName, *inGraphFileName, *outIVfileName ;
double   alpha, nA21, nfent1, nfops1, nL11, nL22, nPhi, nV, t1, t2 ;
Graph    *graph ;
int      ii, inside, J, K, msglvl, nfind1, nfront, nJ, nleaves1, 
         nnode1, nvtx, rc, sizeJ, totalgain, vsize, v, w ;
int      *adjJ, *compids, *nodwghts, *vadj, *vtxToFront, *vwghts ;
IV       *compidsIV ;
IVL      *symbfacIVL ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 7 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inGraphFile alpha"
"\n         outIVfile "
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
"\n    alpha        -- weight parameter"
"\n       alpha = 0 --> minimize storage"
"\n       alpha = 1 --> minimize solve ops"
"\n    outIVfile    -- output file for oldToNew vector,"
"\n                    must be *.ivf or *.ivb"
"\n", argv[0]) ;
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
inETreeFileName  = argv[3] ;
inGraphFileName  = argv[4] ;
alpha            = atof(argv[5]) ;
outIVfileName    = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inGraphFile   -- %s" 
        "\n alpha         -- %f" 
        "\n outIVfile     -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inETreeFileName, inGraphFileName, alpha, outIVfileName) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
if ( strcmp(inETreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
etree = ETree_new() ;
MARKTIME(t1) ;
rc = ETree_readFromFile(etree, inETreeFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in etree from file %s",
        t2 - t1, inETreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_readFromFile(%p,%s)",
           rc, etree, inETreeFileName) ;
   exit(-1) ;
}
ETree_leftJustify(etree) ;
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
nfront     = ETree_nfront(etree) ;
tree       = ETree_tree(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
graph = Graph_new() ;
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
nvtx = graph->nvtx ;
vwghts = graph->vwghts ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Graph object from file %s",
        inGraphFileName) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------
   compute the statistics
   ----------------------
*/
nnode1 = etree->tree->n ;
nfind1 = ETree_nFactorIndices(etree) ;
nfent1 = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
nfops1 = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
nleaves1 = Tree_nleaves(etree->tree) ;
fprintf(stdout, "\n root front %d has %d vertices",
        etree->tree->root,
        etree->nodwghtsIV->vec[etree->tree->root]) ;
/*
   ---------------------------------
   create the symbolic factorization
   ---------------------------------
*/
symbfacIVL = SymbFac_initFromGraph(etree, graph) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
} else {
   IVL_writeStats(symbfacIVL, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------
   find the optimal partition
   --------------------------
*/
compidsIV = ETree_optPart(etree, graph, symbfacIVL, alpha,
                          &totalgain, msglvl, msgFile) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(compidsIV, msgFile) ;
} else {
   IV_writeStats(compidsIV, msgFile) ;
}
fflush(msgFile) ;
compids = IV_entries(compidsIV) ;
/*
   ------------------------------------------------------
   compute the number of vertices in the schur complement
   ------------------------------------------------------
*/
for ( J = 0, nPhi = nV = 0. ; J < nfront ; J++ ) {
   if ( compids[J] == 0 ) {
      nPhi += nodwghts[J] ;
   }
   nV += nodwghts[J] ;
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
      if ( K > J && compids[K] == compids[J] ) {
         inside += (vwghts == NULL) ? 1 : vwghts[w] ;
         if ( msglvl > 3 ) {
            fprintf(msgFile, ", inside") ;
         }
      }
   }
   if ( compids[J] != 0 ) {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    inside = %d, adding %d to L11",
                 inside, nJ*nJ + 2*nJ*inside) ;
      }
      nL11 += (nJ*(nJ+1))/2 + nJ*inside ;
   } else {
      if ( msglvl > 3 ) {
         fprintf(msgFile, "\n    inside = %d, adding %d to L22",
                 inside, (nJ*(nJ+1))/2 + nJ*inside) ;
      }
      nL22 += (nJ*(nJ+1))/2 + nJ*inside ;
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n |L| = %.0f, |L11| = %.0f, |L22| = %.0f",
           nfent1, nL11, nL22) ;
}
/*
   ------------------------------------
   compute the number of entries in A21
   ------------------------------------
*/
nA21 = 0 ;
if ( vwghts != NULL ) {
   for ( v = 0 ; v < nvtx ; v++ ) {
      J = vtxToFront[v] ;
      if ( compids[J] != 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            K = vtxToFront[w] ;
            if ( compids[K] == 0 ) {
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
      J = vtxToFront[v] ;
      if ( compids[J] != 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            K = vtxToFront[w] ;
            if ( compids[K] == 0 ) {
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
           nfent1, nL11, nL22, nA21) ;
   fprintf(msgFile,
      "\n storage: explicit = %.0f, semi-implicit = %.0f, ratio = %.3f"
      "\n opcount: explicit = %.0f, semi-implicit = %.0f, ratio = %.3f",
      nfent1, nL11 + nA21 + nL22,
      nfent1/(nL11 + nA21 + nL22),
      2*nfent1, 4*nL11 + 2*nA21 + 2*nL22,
      2*nfent1/(4*nL11 + 2*nA21 + 2*nL22)) ;
   fprintf(msgFile, "\n ratios %8.3f %8.3f %8.3f",
           nPhi/nV,
           nfent1/(nL11 + nA21 + nL22),
           2*nfent1/(4*nL11 + 2*nA21 + 2*nL22)) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
Graph_free(graph) ;
IVL_free(symbfacIVL) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
