/*  orderViaStages.c  */

#include "../MSMD.h"
#include "../../DSTree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------------
   order a graph using the multi-stage minimum degree algorithm.
   we read in
      (1) a Graph object
      (2) a IV object that contains the stages for the vertices.
      (3) ordering parameters
   and output (optionally)
      (1) and ETree front tree object
      (2) old-to-new IV object
      (3) new-to-old IV object
      
   created -- 96nov01, cca
   -------------------------------------------------------------
*/
{
char        *inGraphFileName, *inIVfileName, *msgFileName, 
            *outETreeFileName, *outNewToOldIVfileName, 
            *outOldToNewIVfileName ;
double      cpu, ops, t1, t2 ;
ETree       *etree ;
FILE        *msgFile ;
Graph       *g ;
int         compressFlag, ierr, msglvl, nfind, nvtx, nzf, 
            prioType, rc, seed, stepType ;
int         *stages ;
IV          *newToOldIV, *oldToNewIV, *stagesIV ;
MSMD        *msmd ;
MSMDinfo    *msmdinfo ;

if ( argc != 12 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inGraphFile inIVfile seed "
"\n        compressFlag prioType stepType oldToNewIVfile "
"\n        newToOldIVfile frontTreeFile " 
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
"\n    inIVfile     -- input file, must be *.ivf or *.ivb"
"\n    seed         -- random number seed"
"\n    compressFlag -- compression flag"
"\n       compressFlag / 4 >= 1 --> compress before elimination"
"\n       compressFlag %% 4 == 2 --> compress at each elimination step,"
"\n                                 consider all nodes"
"\n       compressFlag %% 4 == 1 --> compress at each elimination step,"
"\n                                 but only consider 2-adj nodes"
"\n       compressFlag %% 4 == 0 --> do not perform any compression"
"\n   prioType -- update type"
"\n      prioType == 1 --> true updates"
"\n      prioType == 2 --> approximate updates"
"\n      prioType == 3 --> half and half"
"\n      prioType == 4 --> maximal independent set"
"\n   stepType -- degree extent for independent set elimination"
"\n      when stepType < 1 --> only one node is eliminated at a step,"
"\n         e.g., like QMD from SPARSPAK and YSMP"
"\n      stepType == 1 --> regular multiple elimination, e.g., GENMMD"
"\n      stepType >  1 --> extended multiple elimination"
"\n       an independent set of nodes is selected for elimination"
"\n       whose degree satisfies mindeg <= degree <= stepType*mindegree"
"\n   oldToNewIVfile -- file for old-to-new IV object"
"\n   newToOldIVfile -- file for new-to-old IV object"
"\n   frontTreeFile  -- file for front tree's ETree object"
"\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
msgFileName = argv[2] ;
if ( strcmp(msgFileName, "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(msgFileName, "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], msgFileName) ;
   return(-1) ;
}
inGraphFileName       = argv[3] ;
inIVfileName          = argv[4] ;
seed                  = atoi(argv[5]) ;
compressFlag          = atoi(argv[6]) ;
prioType              = atoi(argv[7]) ;
stepType              = atoi(argv[8]) ;
outETreeFileName      = argv[9] ;
outOldToNewIVfileName = argv[10] ;
outNewToOldIVfileName = argv[11] ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl          -- %d" 
        "\n msgFile         -- %s" 
        "\n inGraphFile     -- %s" 
        "\n inIVfile        -- %s" 
        "\n seed            -- %d" 
        "\n compressFlag    -- %d" 
        "\n prioType        -- %d" 
        "\n stepType        -- %d" 
        "\n outETreeFile    -- %s" 
        "\n outOldToNewFile -- %s" 
        "\n outNewToOldFile -- %s" 
        "\n", argv[0], msglvl, msgFileName, inGraphFileName, 
        inIVfileName, seed, 
        compressFlag, prioType, stepType, outETreeFileName,
        outOldToNewIVfileName, outNewToOldIVfileName) ;
fflush(msgFile) ;
/*
   ------------------------------------------
   initialize the MSMDinfo information object
   ------------------------------------------
*/
msmdinfo                = MSMDinfo_new() ;
msmdinfo->seed          = seed ;
msmdinfo->compressFlag  = compressFlag ;
msmdinfo->prioType      = prioType ;
msmdinfo->stepType      = stepType ;
msmdinfo->msglvl        = msglvl        ;
msmdinfo->msgFile       = msgFile       ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
g = Graph_new() ;
Graph_setDefaultFields(g) ;
if ( (rc = Graph_readFromFile(g, inGraphFileName)) != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
        rc, g, inGraphFileName) ;
   exit(-1) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s", 
        t2 - t1, inGraphFileName) ;
nvtx = g->nvtx ;
if ( msglvl < 4 ) {
   Graph_writeStats(g, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeForHumanEye(g, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the IV object
   -------------------------
*/
stagesIV = IV_new() ;
IV_setDefaultFields(stagesIV) ;
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
} else {
   MARKTIME(t1) ;
   if ( (rc = IV_readFromFile(stagesIV, inIVfileName)) != 1 ) {
      fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
              rc, stagesIV, inIVfileName) ;
      exit(-1) ;
   }
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : read in IV object from file %s", 
           t2 - t1, inIVfileName) ;
}
if ( msglvl < 4 ) {
   IV_writeStats(stagesIV, msgFile) ;
   fflush(msgFile) ;
} else {
   IV_writeForHumanEye(stagesIV, msgFile) ;
   fflush(msgFile) ;
}
stages = IV_entries(stagesIV) ;
/*
   --------------
   order via msmd
   --------------
*/
MARKTIME(t1) ;
msmd = MSMD_new() ;
MSMD_order(msmd, g, stages, msmdinfo) ;
MARKTIME(t2) ;
cpu = t2 - t1 ;
fprintf(msgFile, "\n CPU %9.5f : order the graph", cpu) ;
fflush(msgFile) ;
MSMDinfo_print(msmdinfo, msgFile) ;
fflush(msgFile) ;
/*
   ----------------------
   extract the front tree
   ----------------------
*/
MARKTIME(t1) ;
etree = MSMD_frontETree(msmd) ;
nfind = ETree_nFactorIndices(etree) ;
nzf   = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
ops   = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
MARKTIME(t2) ;
fprintf(msgFile,  
        "\n\n CPU %9.5f : make the front tree", t2 - t1) ;
fprintf(msgFile,  
        "\n FACTOR : %9d indices, %9d entries, %9.0f operations", 
        nfind, nzf, ops) ;
if ( msglvl < 3 ) {
   ETree_writeStats(etree, msgFile) ;
   fflush(msgFile) ;
} else {
   ETree_writeForHumanEye(etree, msgFile) ;
   fflush(msgFile) ;
}
fprintf(msgFile, "\n STATS    %10d %10.0f %8.3f", nzf, ops, cpu) ;
/*
   -------------------------------------------
   write the front tree to a file if requested
   -------------------------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   ETree_writeToFile(etree, outETreeFileName) ;
}
/*
   -----------------------------------------------------------------
   generate the permutation vectors and write to a file if requested
   -----------------------------------------------------------------
*/
if ( strcmp(outOldToNewIVfileName, "none") != 0 ) {
   oldToNewIV = IV_new() ;
} else {
   oldToNewIV = NULL ;
}
if ( strcmp(outNewToOldIVfileName, "none") != 0 ) {
   newToOldIV = IV_new() ;
} else {
   newToOldIV = NULL ;
}
if ( oldToNewIV != NULL || newToOldIV != NULL ) {
   MSMD_fillPerms(msmd, newToOldIV, oldToNewIV) ;
}
if ( oldToNewIV != NULL ) {
   IV_writeToFile(oldToNewIV, outOldToNewIVfileName) ;
}
if ( newToOldIV != NULL ) {
   IV_writeToFile(newToOldIV, outNewToOldIVfileName) ;
}
/*
   ----------------------------
   free all the working storage
   ----------------------------
*/
Graph_free(g) ;
MSMD_free(msmd) ;
MSMDinfo_free(msmdinfo) ;
ETree_free(etree) ;
if ( newToOldIV != NULL ) {
   IV_free(newToOldIV) ;
}
if ( oldToNewIV != NULL ) {
   IV_free(oldToNewIV) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(-1) ; }

/*--------------------------------------------------------------------*/
