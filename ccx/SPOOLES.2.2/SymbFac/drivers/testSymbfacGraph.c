/*  testSymbfacGraph.c  */

#include "../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   (1) read in an ETree object.
   (2) get the old-to-new vertex permutation.
   (3) permute the vtx-to-front map.
   (4) get the symbolic factorization IVL object.
   (5) permute the ETree object.
   (6) optionally write the permuted ETree object to a file
   (7) optionally write the old-to-new IV object to a file
   (8) optionally write the symbolic factorization IV object to a file
   
   created -- 96oct03, cca
   ---------------------------------------------------------------
*/
{
char     *inETreeFileName, *inGraphFileName, *outETreeFileName,
         *outIVfileName, *outIVLfileName ;
double   nfops1, t1, t2 ;
Graph    *graph ;
int      msglvl, nfent1, nfind1, nleaves1, nnode1, rc ;
IV       *vtxOldToNewIV ;
IVL      *symbfacIVL ;
ETree    *etree ;
FILE     *msgFile ;

if ( argc != 8 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inGraphFile outETreeFile"
"\n         outIVfile outIVLfile"
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
"\n    outETreeFile -- output file, must be *.etreef or *.etreeb"
"\n    outIVfile    -- output file for oldToNew vector,"
"\n                    must be *.ivf or *.ivb"
"\n    outIVLfile   -- output file for symbolic factorization object"
"\n                    must be *.ivlf or *.ivlb"
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
outETreeFileName = argv[5] ;
outIVfileName    = argv[6] ;
outIVLfileName   = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inGraphFile   -- %s" 
        "\n outETreeFile  -- %s" 
        "\n outIVfile     -- %s" 
        "\n outIVLfile    -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inETreeFileName, inGraphFileName, outETreeFileName, 
        outIVfileName, outIVLfileName) ;
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
/*
ETree_leftJustify(etree) ;
*/
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
/*
   int   front, nfront, nvtx, v ;
   int   *head, *link, *vtxToFront ;

   nfront = etree->nfront ;
   nvtx   = etree->nvtx   ;
   head = IVinit(nfront, -1) ;
   link = IVinit(nvtx, -1) ;
   vtxToFront = IV_entries(etree->vtxToFrontIV) ;
   for ( v = nvtx - 1 ; v >= 0 ; v-- ) {
      front = vtxToFront[v] ;
      link[v] = head[front] ;
      head[front] = v ;
   }
   for ( front = 0 ; front < nfront ; front++ ) {
      fprintf(msgFile, "\n front %3d :", front) ;
      for ( v = head[front] ; v != -1 ; v = link[v] ) {
         fprintf(msgFile, " %d", v) ;
      }
   }
   IVfree(head) ;
   IVfree(link) ;
*/
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
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
   -----------------------------
   get the permutation IV object
   -----------------------------
*/
vtxOldToNewIV = ETree_oldToNewVtxPerm(etree) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   IV_writeForHumanEye(vtxOldToNewIV, msgFile) ;
} else {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   IV_writeStats(vtxOldToNewIV, msgFile) ;
}
fflush(msgFile) ;
IV_writeToFile(vtxOldToNewIV, "oldToNew.ivf") ;
/*
   ----------------------------------
   optionally write out the IV object
   ----------------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(vtxOldToNewIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write vtxOldToNewIV to file %s",
           t2 - t1, outIVfileName) ;
}
/*
   --------------------------------------------
   create the symbolic factorization IVL object
   --------------------------------------------
*/
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromGraph(etree, graph) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compute the symbolic factorization",
        t2 - t1) ;
fprintf(msgFile, 
        "\n\n symbolic factorization IVL object in old ordering") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
} else {
   IVL_writeStats(symbfacIVL, msgFile) ;
}
fflush(msgFile) ;
MARKTIME(t1) ;
IVL_overwrite(symbfacIVL, vtxOldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : permute symbfac", t2 - t1) ;
fprintf(msgFile, 
        "\n\n symbolic factorization IVL object after overwrite") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
} else {
   IVL_writeStats(symbfacIVL, msgFile) ;
}
fflush(msgFile) ;
MARKTIME(t1) ;
IVL_sortUp(symbfacIVL) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : sort up", t2 - t1) ;
fprintf(msgFile, 
        "\n\n symbolic factorization IVL object in new ordering") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
} else {
   IVL_writeStats(symbfacIVL, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------------------------
   permute the vertices in the ETree object
   ----------------------------------------
*/
fprintf(msgFile, "\n\n before permuting the vertices") ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
MARKTIME(t1) ;
ETree_permuteVertices(etree, vtxOldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : permute vertices in ETree", t2 - t1) ;
fprintf(msgFile, "\n\n after permuting the vertices") ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------------------
   optionally write out the ETree object
   -------------------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, etree, outETreeFileName) ;
}
/*
   -----------------------------------
   optionally write out the IVL object
   -----------------------------------
*/
if ( strcmp(outIVLfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IVL_writeToFile(symbfacIVL, outIVLfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write symbfac IVL to file %s",
           t2 - t1, outIVLfileName) ;
}
/*
{
int   count, ii, J, nfront, nvtx, sizeJ, v, w ;
int   *head, *indJ, *link, *vtxToFront ;

   nvtx   = graph->nvtx ;
   nfront = etree->nfront ;
   head   = IVinit(nfront, -1) ;
   link   = IVinit(nvtx, -1) ;
   vtxToFront = ETree_vtxToFront(etree) ;
   for ( v = nvtx - 1 ; v >= 0 ; v-- ) {
      J = vtxToFront[v] ;
      link[v] = head[J] ;
      head[J] = v ;
   }
   fprintf(msgFile, "\n /adjncy [") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      IVL_listAndSize(symbfacIVL, J, &sizeJ, &indJ) ;
      for ( v = head[J] ; v != -1 ; v = link[v] ) {
         fprintf(msgFile, "\n [") ;
         for ( ii = 0 ; ii < sizeJ ; ii++ ) {
            w = indJ[ii] ;
            if ( v <= w ) {
               fprintf(msgFile, " %d", w) ;
            }
         }
         fprintf(msgFile, " ]") ;
      }
   }
   fprintf(msgFile, "\n ] def") ;
   fprintf(msgFile, "\n /fsinfo [") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      fprintf(msgFile, "\n [ %d", head[J]) ;
      for ( v = head[J], count = 0 ; v != -1 ; v = link[v] ) {
         count++ ;
      }
      fprintf(msgFile, " %d ]", count) ;
   }
   fprintf(msgFile, "\n ] def") ;
}
*/
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
Graph_free(graph) ;
IV_free(vtxOldToNewIV) ;
IVL_free(symbfacIVL) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
