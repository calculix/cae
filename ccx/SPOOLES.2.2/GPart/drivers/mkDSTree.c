/*  mkDSTree.c  */

#include "../GPart.h"
#include "../../DSTree.h"
#include "../../MSMD.h"
#include "../../BKL.h"
#include "../../ETree.h"
#include "../../Perm.h"
#include "../../IV.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   test the recursive bisection algorithm that uses
   (1) fishnet to get the domain decomposition
   (2) domain/segment BKL to get the two set partition
   (3) Dulmadge-Mendelsohn decomposition to smooth the bisector
   the output is a DSTree object to hold the domain/separator tree

   created -- 96mar09, cca
   ---------------------------------------------------------------
*/
{
char        *inGraphFileName, *msgFileName, *outDSTreeFileName ;
DSTree      *dstree ;
DDsepInfo   *info ;
double      alpha, freeze, msCPU, t1, t2 ;
FILE        *msgFile ;
GPart       *gpart ;
Graph       *gf ;
int         DDoption, maxdomweight, maxweight, minweight, msglvl, 
            nlayer, nvtx, rc, seed ;

if ( argc != 13 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile seed"
      "\n         minweight maxweight freeze alpha maxdomwght "
      "\n         DDoption nlayer outDSTreeFileName"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
      "\n    seed         -- random number seed"
      "\n    minweight    -- minimum domain weight"
      "\n    maxweight    -- maximum domain weight"
      "\n    freeze       -- cutoff multiplier for nodes of high degree"
      "\n    alpha        -- cost function parameter"
      "\n    maxdomweight -- maximum subgraph weight"
      "\n    DDoption     -- option for domain decomposition"
      "\n       1 --> fishnet for each subgraph"
      "\n       2 --> fishnet for graph, projection for each subgraph"
      "\n    nlayer -- number of layers for max flow improvement"
      "\n    outDSTreeFileName -- output file, must be *.dstreef"
      "\n                         or *.dstreeb"
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
inGraphFileName   = argv[3] ;
seed              = atoi(argv[4]) ;
minweight         = atoi(argv[5]) ;
maxweight         = atoi(argv[6]) ;
freeze            = atof(argv[7]) ;
alpha             = atof(argv[8]) ;
maxdomweight      = atoi(argv[9]) ;
DDoption          = atoi(argv[10]) ;
nlayer            = atoi(argv[11]) ;
outDSTreeFileName = argv[12] ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inGraphFile   -- %s" 
        "\n seed          -- %d" 
        "\n minweight     -- %d" 
        "\n maxweight     -- %d" 
        "\n freeze        -- %f" 
        "\n alpha         -- %f" 
        "\n maxdomweight  -- %d" 
        "\n DDoption      -- %d" 
        "\n nlayer        -- %d" 
        "\n outDSTreeFile -- %s" 
        "\n", argv[0], msglvl, msgFileName, inGraphFileName, seed, 
       minweight, maxweight, freeze, alpha, maxdomweight, DDoption,
       nlayer, outDSTreeFileName) ;
fflush(msgFile) ;
/*
   ---------------------------------------
   initialize the DDsep information object
   ---------------------------------------
*/
info                = DDsepInfo_new() ;
info->seed          = seed ;
info->minweight     = minweight ;
info->maxweight     = maxweight ;
info->freeze        = freeze ;
info->alpha         = alpha ;
info->DDoption      = DDoption ;
info->maxcompweight = maxdomweight ;
info->nlayer        = nlayer        ;
info->msglvl        = msglvl        ;
info->msgFile       = msgFile       ;
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
gf = Graph_new() ;
Graph_setDefaultFields(gf) ;
if ( (rc = Graph_readFromFile(gf, inGraphFileName)) != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
        rc, gf, inGraphFileName) ;
   exit(-1) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s", 
        t2 - t1, inGraphFileName) ;
nvtx = gf->nvtx ;
if ( msglvl < 4 ) {
   Graph_writeStats(gf, msgFile) ;
   fflush(msgFile) ;
} else {
   Graph_writeForHumanEye(gf, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------
   create the GPart object
   -----------------------
*/
MARKTIME(t1) ;
gpart = GPart_new() ;
GPart_init(gpart, gf) ;
GPart_setMessageInfo(gpart, msglvl, msgFile) ;
MARKTIME(t2) ;
/*
   ------------------------------------------
   get the DSTree object that represents the
   domain/separator partition of the vertices
   ------------------------------------------
*/
MARKTIME(t1) ;
dstree = GPart_RBviaDDsep(gpart, info) ;
MARKTIME(t2) ;
msCPU = t2 - t1 ;
fprintf(msgFile, "\n\n CPU %9.5f : find subgraph tree ", msCPU) ;
DDsepInfo_writeCpuTimes(info, msgFile) ;

if ( msglvl > 0 ) {
   fprintf(msgFile, "\n # subgraphs = %d", dstree->tree->n) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n DSTree subgraph tree") ;
   DSTree_writeForHumanEye(dstree, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n map from vertices to subgraphs") ;
   IV_writeForHumanEye(dstree->mapIV, msgFile) ;
   fflush(msgFile) ;
}
DDsepInfo_free(info) ;
/*
   --------------------------------------------
   renumber the tree via a post-order traversal
   --------------------------------------------
*/
DSTree_renumberViaPostOT(dstree) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n renumbered DSTree subgraph tree") ;
   DSTree_writeForHumanEye(dstree, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------
   optionally write the DSTree object to a file
   --------------------------------------------
*/
if ( strcmp(outDSTreeFileName, "none") != 0 ) {
   DSTree_writeToFile(dstree, outDSTreeFileName) ;
}
/*
   ----------------------------
   free all the working storage
   ----------------------------
*/
Graph_free(gpart->g) ;
GPart_free(gpart) ;
DSTree_free(dstree) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
