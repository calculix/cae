/*  testRBviaDDsep2.c  */

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
   ------------------------------------------------------------
   test the recursive bisection algorithm that uses
   (1) fishnet to get the domain decomposition
   (2) domain/segment BKL to get the two set partition
   (3) Dulmadge-Mendelsohn decomposition to smooth the bisector

   created -- 96mar09, cca
   ------------------------------------------------------------
*/
{
char        *inGraphFileName, *msgFileName ;
DSTree      *dstree ;
DDsepInfo   *info ;
double      alpha, freeze, msCPU, msops, ndCPU, ndops, 
            phiFrac, rbCPU, t1, t2 ;
double      *fracvec, *mscpuvec, *msopsvec, 
            *ndcpuvec, *ndopsvec, *rbcpuvec ;
ETree       *etree ;
FILE        *msgFile ;
GPart       *gpart ;
Graph       *gf ;
int         DDoption, maxdomweight, maxweight, minweight, 
            msnfind, msnzf, msglvl, 
            ndnfind, ndnzf, nlayer, ntest, nvtx, phiWeight, rc, seed ;
int         *emap, *msnzfvec, *ndnzfvec ;
IV          *stagesIV ;
MSMD        *msmd ;
MSMDinfo    *msmdinfo ;

if ( argc != 12 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inGraphFile ntest"
      "\n         minweight maxweight freeze alpha maxdomwght "
      "\n         DDoption nlayer"
      "\n    msglvl       -- message level"
      "\n    msgFile      -- message file"
      "\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
      "\n    ntest        -- number of tests to make"
      "\n    minweight    -- minimum domain weight"
      "\n    maxweight    -- maximum domain weight"
      "\n    freeze       -- cutoff multiplier for nodes of high degree"
      "\n    alpha        -- cost function parameter"
      "\n    maxdomweight -- maximum subgraph weight"
      "\n    DDoption     -- option for domain decomposition"
      "\n       1 --> fishnet for each subgraph"
      "\n       2 --> fishnet for graph, projection for each subgraph"
      "\n    nlayer -- number of layers for max flow improvement"
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
inGraphFileName = argv[3] ;
ntest           = atoi(argv[4]) ;
if ( ntest % 2 == 0 ) {
   ntest++ ;
}
minweight       = atoi(argv[5]) ;
maxweight       = atoi(argv[6]) ;
freeze          = atof(argv[7]) ;
alpha           = atof(argv[8]) ;
maxdomweight    = atoi(argv[9]) ;
DDoption        = atoi(argv[10]) ;
nlayer          = atoi(argv[11]) ;
fprintf(msgFile, 
        "\n %s : "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inGraphFile  -- %s" 
        "\n ntest        -- %d" 
        "\n minweight    -- %d" 
        "\n maxweight    -- %d" 
        "\n freeze       -- %f" 
        "\n alpha        -- %f" 
        "\n maxdomweight -- %d" 
        "\n DDoption     -- %d" 
        "\n nlayer       -- %d" 
        "\n", argv[0], msglvl, msgFileName, inGraphFileName, ntest, 
       minweight, maxweight, freeze, alpha, maxdomweight, DDoption,

       nlayer) ;
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
if ( msglvl > 3 ) {
   Graph_writeForHumanEye(gf, msgFile) ;
   fflush(msgFile) ;
} else if ( msglvl > 1 ) {
   Graph_writeStats(gf, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------
   allocate the statistics vectors
   -------------------------------
*/
fracvec  = DVinit(ntest, 0.0) ;
ndnzfvec = IVinit(ntest,  0)  ;
msnzfvec = IVinit(ntest,  0)  ;
ndopsvec = DVinit(ntest, 0.0) ;
msopsvec = DVinit(ntest, 0.0) ;
rbcpuvec = DVinit(ntest, 0.0) ;
ndcpuvec = DVinit(ntest, 0.0) ;
mscpuvec = DVinit(ntest, 0.0) ;
/*
   -------------------
   loop over the tests
   -------------------
*/
for ( seed = 1 ; seed <= ntest ; seed++ ) {
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
   rbCPU = t2 - t1 ;
   phiWeight = DSTree_separatorWeight(dstree, gf->vwghts) ;
   fprintf(msgFile, 
       "\n\n CPU %9.5f : find subgraph tree, %d subgraphs, |Phi| = %d ",
       rbCPU, dstree->tree->n, phiWeight) ;
   if ( msglvl > 0 ) {
       DDsepInfo_writeCpuTimes(info, msgFile) ;
   }
/*
   --------------------------------------------
   compute the weight of the separator vertices
   --------------------------------------------
*/
   if ( gf->totvwght == 0 ) {
      if ( gf->vwghts == NULL ) {
         gf->totvwght = nvtx ;
      } else {
         gf->totvwght = IVsum(nvtx, gf->vwghts) ;
      }
   }
   phiFrac = ((double) phiWeight) / gf->totvwght ;
   
   if ( msglvl > 1 ) {
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
   -----------------------------------------
   set the stages for two-stage multisection
   -----------------------------------------
*/
   stagesIV = DSTree_MS2stages(dstree) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n stages for MS2") ;
      IV_writeForHumanEye(stagesIV, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   read in the graph 
   -----------------
   Graph_readFromFile(gf, inGraphFileName) ;
*/
/*
   -------------------------------
   order as two-stage multisection
   -------------------------------
*/
   msmdinfo = MSMDinfo_new() ;
   msmdinfo->compressFlag = 2       ;
   msmdinfo->prioType     = 1       ;
   msmdinfo->stepType     = 1       ;
   msmdinfo->seed         = seed    ;
   msmdinfo->msglvl       = msglvl  ;
   msmdinfo->msgFile      = msgFile ;
   MARKTIME(t1) ;
   msmd = MSMD_new() ;
   MSMD_order(msmd, gf, IV_entries(stagesIV), msmdinfo) ;
   MARKTIME(t2) ;
   msCPU = t2 - t1 ;
   fprintf(msgFile, "\n CPU %9.5f : order the graph via MS", msCPU) ;
   fflush(msgFile) ;
   if ( msglvl > 1 ) {
      MSMDinfo_print(msmdinfo, msgFile) ;
      fflush(msgFile) ;
   }
   IV_free(stagesIV) ;
/*
   ----------------------
   extract the front tree
   ----------------------
*/
   MARKTIME(t1) ;
   emap    = IVinit(nvtx, -1) ;
   etree   = MSMD_frontETree(msmd) ;
   msnfind = ETree_nFactorIndices(etree) ;
   msnzf   = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
   msops   = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   MARKTIME(t2) ;
   fprintf(msgFile,  
           "\n CPU %9.5f : make the elimination tree", t2 - t1) ;
   fprintf(msgFile,  
           "\n MS FACTOR : %9d indices, %9d entries, %9.0f operations", 
           msnfind, msnzf, msops) ;
   if ( msglvl > 3 ) {
      ETree_writeForHumanEye(etree, msgFile) ;
      fflush(msgFile) ;
   } else if ( msglvl > 1 ) {
      ETree_writeStats(etree, msgFile) ;
      fflush(msgFile) ;
   }
   fprintf(msgFile, "\n STATSMS2 %10d %10.0f %8.3f %8.3f %8.3f",
           msnzf, msops, rbCPU, msCPU, rbCPU + msCPU) ;
   MSMD_free(msmd) ;
   MSMDinfo_free(msmdinfo) ;
   ETree_free(etree) ;
   IVfree(emap) ;
/*
   ------------------------------------
   set the stages for nested dissection
   ------------------------------------
*/
   stagesIV = DSTree_NDstages(dstree) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n stages for ND") ;
      IV_writeForHumanEye(stagesIV, msgFile) ;
      fflush(msgFile) ;
   }
/*
   -----------------
   read in the graph 
   -----------------
   Graph_readFromFile(gf, inGraphFileName) ;
*/
/*
   -------------------------------------
   order as incomplete nested dissection
   -------------------------------------
*/
   msmdinfo = MSMDinfo_new() ;
   msmdinfo->compressFlag = 2       ;
   msmdinfo->prioType     = 1       ;
   msmdinfo->stepType     = 1       ;
   msmdinfo->seed         = seed    ;
   msmdinfo->msglvl       = msglvl  ;
   msmdinfo->msgFile      = msgFile ;
   MARKTIME(t1) ;
   msmd = MSMD_new() ;
   MSMD_order(msmd, gf, IV_entries(stagesIV), msmdinfo) ;
   MARKTIME(t2) ;
   ndCPU = t2 - t1 ;
   fprintf(msgFile, "\n CPU %9.5f : order the graph via ND", ndCPU) ;
   fflush(msgFile) ;
   if ( msglvl > 1 ) {
      MSMDinfo_print(msmdinfo, msgFile) ;
      fflush(msgFile) ;
   }
   IV_free(stagesIV) ;
/*
   ----------------------
   extract the front tree
   ----------------------
*/
   MARKTIME(t1) ;
   emap    = IVinit(nvtx, -1) ;
   etree   = MSMD_frontETree(msmd) ;
   ndnfind = ETree_nFactorIndices(etree) ;
   ndnzf   = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
   ndops   = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   MARKTIME(t2) ;
   fprintf(msgFile,  
           "\n CPU %9.5f : make the elimination tree", t2 - t1) ;
   fprintf(msgFile,  
           "\n ND FACTOR : %9d indices, %9d entries, %9.0f operations", 
           ndnfind, ndnzf, ndops) ;
   if ( msglvl > 3 ) {
      ETree_writeForHumanEye(etree, msgFile) ;
      fflush(msgFile) ;
   } else if ( msglvl > 1 ) {
      ETree_writeStats(etree, msgFile) ;
      fflush(msgFile) ;
   }
   fprintf(msgFile, "\n STATSND  %10d %10.0f %8.3f %8.3f %8.3f",
           ndnzf, ndops, rbCPU, ndCPU, rbCPU + ndCPU) ;
   MSMD_free(msmd) ;
   MSMDinfo_free(msmdinfo) ;
   ETree_free(etree) ;
   IVfree(emap) ;
/*
   ------------------------
   print out the statistics
   ------------------------
*/
   fprintf(msgFile, 
           "\n ALL %6.3f %8.3f %8d %10.0f %8.3f %8d %10.0f %8.3f",
           phiFrac, rbCPU, ndnzf, ndops, ndCPU, msnzf, msops, msCPU) ;
/*
   ----------------------------
   free all the working storage
   ----------------------------
*/
   GPart_free(gpart) ;
   DSTree_free(dstree) ;
/*
   --------------------
   store the statistics
   --------------------
*/
   fracvec[seed-1]  = phiFrac ;
   rbcpuvec[seed-1] = rbCPU   ;
   ndcpuvec[seed-1] = ndCPU   ;
   mscpuvec[seed-1] = msCPU   ;
   ndnzfvec[seed-1] = ndnzf   ;
   msnzfvec[seed-1] = msnzf   ;
   ndopsvec[seed-1] = ndops   ;
   msopsvec[seed-1] = msops   ;
}
/*
   ---------------------------
   print out the median values
   ---------------------------
*/
DVqsortUp(ntest, fracvec) ;
DVqsortUp(ntest, rbcpuvec) ;
DVqsortUp(ntest, ndcpuvec) ;
DVqsortUp(ntest, mscpuvec) ;
IVqsortUp(ntest, ndnzfvec) ;
IVqsortUp(ntest, msnzfvec) ;
DVqsortUp(ntest, ndopsvec) ;
DVqsortUp(ntest, msopsvec) ;
fprintf(msgFile, 
        "\n MEDIAN %6.3f %8.3f %8d %10.0f %8.3f %8d %10.0f %8.3f",
        fracvec[ntest/2], rbcpuvec[ntest/2], ndnzfvec[ntest/2], 
        ndopsvec[ntest/2],   ndcpuvec[ntest/2], msnzfvec[ntest/2], 
        msopsvec[ntest/2],   mscpuvec[ntest/2]) ;
/*
   ----------------------------
   free all the working storage
   ----------------------------
*/
Graph_free(gf) ;
DVfree(fracvec) ;
IVfree(ndnzfvec) ;
IVfree(msnzfvec) ;
DVfree(ndopsvec) ;
DVfree(msopsvec) ;
DVfree(rbcpuvec) ;
DVfree(ndcpuvec) ;
DVfree(mscpuvec) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
