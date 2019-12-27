/*  createETree.c  */

#include "../ETree.h"
#include "../../Perm.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------
   read in a Graph and a Perm object.
   create the ETree object and fill an IV object 
   with the compids of the two-set partition

   created -- 96may02, cca
   ---------------------------------------------
*/
{
char     *inGraphFileName, *inPermFileName, 
         *outETreeFileName, *outIVfileName ;
double   t1, t2 ;
int      msglvl, rc ;
ETree    *etree, *fsETree ;
IV       *fsMapIV ;
FILE     *msgFile ;
Graph    *graph ;
Perm     *perm ;

if ( argc != 7 ) {
   fprintf(stdout, 
   "\n\n usage : %s msglvl msgFile inGraphFile inPermFile "
   "\n         outIVfile outETreeFile"
   "\n    msglvl       -- message level"
   "\n    msgFile      -- message file"
   "\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
   "\n    inPermFile   -- input file, must be *.permf or *.permb"
   "\n    outIVfile    -- output file for compids[]"
   "\n                    must be *.ivf or *.ivf"
   "\n    outETreeFile -- output file, must be *.etreef or *.etreeb"
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
inGraphFileName  = argv[3] ;
inPermFileName   = argv[4] ;
outIVfileName    = argv[5] ;
outETreeFileName = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inGraphFile  -- %s" 
        "\n inPermFile   -- %s" 
        "\n outIVfile    -- %s" 
        "\n outETreeFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inGraphFileName, inPermFileName, 
        outIVfileName, outETreeFileName) ;
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
   ------------------------
   read in the Perm object
   ------------------------
*/
if ( strcmp(inPermFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
/*
   exit(0) ;
*/
   perm = NULL ;
} else {
   perm = Perm_new() ;
   MARKTIME(t1) ;
   rc = Perm_readFromFile(perm, inPermFileName) ;
   Perm_fillOldToNew(perm) ;
   Perm_fillNewToOld(perm) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : read in perm from file %s",
           t2 - t1, inPermFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from Perm_readFromFile(%p,%s)",
              rc, perm, inPermFileName) ;
      exit(-1) ;
   }
   rc = Perm_checkPerm(perm) ;
   if ( rc != 1 ) {
      fprintf(stderr, "\n fatal error, Perm not valid") ;
      Perm_writeForHumanEye(perm, stderr) ;
      exit(0) ;
   }
   fprintf(msgFile, "\n\n after reading Perm object from file %s",
           inPermFileName) ;
   if ( msglvl > 0 ) {
      Perm_writeForHumanEye(perm, msgFile) ;
   } else {
      Perm_writeStats(perm, msgFile) ;
   }
   fflush(msgFile) ;
}
fprintf(msgFile, "\n newToOld") ;
IVfprintf(msgFile, perm->size, perm->newToOld) ;
fprintf(msgFile, "\n oldToNew") ;
IVfprintf(msgFile, perm->size, perm->oldToNew) ;
/*
   -----------------------
   create the ETree object
   -----------------------
*/
etree = ETree_new() ;
if ( perm == NULL ) {
   ETree_initFromGraph(etree, graph) ;
} else {
   ETree_initFromGraphWithPerms(etree, graph, perm->newToOld,
                                perm->oldToNew) ;
}
fprintf(msgFile, "\n\n vertex etree") ;
fprintf(msgFile, "\n %d factor indices",
        ETree_nFactorIndices(etree)) ;
fprintf(msgFile, "\n symmetric: %d factor entries",
        ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n nonsymmetric: %d factor entries",
        ETree_nFactorEntries(etree, SPOOLES_NONSYMMETRIC)) ;
fprintf(msgFile, "\n real symmetric       : %.0f factor operations",
        ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n real nonsymmetric    : %.0f factor operations",
        ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_NONSYMMETRIC)) ;
fprintf(msgFile, "\n complex symmetric    : %.0f factor operations",
        ETree_nFactorOps(etree, SPOOLES_COMPLEX, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n complex nonsymmetric : %.0f factor operations",
        ETree_nFactorOps(etree, SPOOLES_COMPLEX, SPOOLES_NONSYMMETRIC));
fsMapIV = ETree_fundSupernodeMap(etree) ;
fsETree = ETree_compress(etree, fsMapIV) ;
fprintf(msgFile, "\n\n fundamental supernode etree") ;
fprintf(msgFile, "\n %d factor indices",
        ETree_nFactorIndices(fsETree)) ;
fprintf(msgFile, "\n symmetric: %d factor entries",
        ETree_nFactorEntries(fsETree, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n nonsymmetric: %d factor entries",
        ETree_nFactorEntries(fsETree, SPOOLES_NONSYMMETRIC)) ;
fprintf(msgFile, "\n real symmetric       : %.0f factor operations",
        ETree_nFactorOps(fsETree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n real nonsymmetric    : %.0f factor operations",
        ETree_nFactorOps(fsETree, SPOOLES_REAL, SPOOLES_NONSYMMETRIC)) ;
fprintf(msgFile, "\n complex symmetric    : %.0f factor operations",
        ETree_nFactorOps(fsETree, SPOOLES_COMPLEX, SPOOLES_SYMMETRIC)) ;
fprintf(msgFile, "\n complex nonsymmetric : %.0f factor operations",
      ETree_nFactorOps(fsETree, SPOOLES_COMPLEX, SPOOLES_NONSYMMETRIC));
fprintf(msgFile, "\n %.0f factor operations",
        ETree_nFactorOps(fsETree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(fsETree, msgFile) ;
} else {
   ETree_writeStats(fsETree, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(fsETree, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, fsETree, outETreeFileName) ;
}
/*
   -------------------------------
   write out the compids IV object
   -------------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(fsETree->vtxToFrontIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outIVfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
           rc, fsETree->vtxToFrontIV, outIVfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
Graph_free(graph)   ;
Perm_free(perm)     ;
ETree_free(etree)   ;
ETree_free(fsETree) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
