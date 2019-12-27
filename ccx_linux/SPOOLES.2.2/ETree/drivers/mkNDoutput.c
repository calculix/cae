/*  mkNDoutput.c  */

#include "../../ETree.h"
#include "../../SymbFac.h"
#include "../../EGraph.h"
#include "../../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   create
   (1) an ETree object for nested dissection on a regular grid
       using a bound on zeros in a front and a bound on front size
   (2) an IV object that maps fronts to threads using a wrap map,
       a balanced map, a subtree-subset map, or a domain
       decomposition map

   created -- 98feb05, cca
   ---------------------------------------------------------------
*/
{
char     *outETreeFileName, *outMapIVfileName ;
double   cutoff, t1, t2 ;
double   ops[5] ;
DV       *cumopsDV ;
int      maptype, maxsize, maxzeros, msglvl, n1, n2, n3, nthread, nvtx, 
         rc, v ;
int      nfind[5], nfronts[5], nzf[5] ;
int      *newToOld, *oldToNew ;
IV       *msIV, *nzerosIV, *ownersIV ;
IVL      *symbfacIVL ;
EGraph   *egraph ;
ETree    *etree0, *etree1, *etree2, *etree3, *etree4 ;
FILE     *msgFile ;
Graph    *graph ;

if ( argc != 13 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile n1 n2 n3 maxzeros maxsize "
      "\n         nthread maptype cutoff etreeFile mapFile"
      "\n    msglvl    -- message level"
      "\n    msgFile   -- message file"
      "\n    n1        -- number of points in the first direction"
      "\n    n2        -- number of points in the second direction"
      "\n    n3        -- number of points in the third direction"
      "\n    maxzeros  -- number of points in the third direction"
      "\n    maxsize   -- maximum number of vertices in a front"
      "\n    nthread   -- number of threads"
      "\n    maptype   -- map type"
      "\n       1 -- wrap map"
      "\n       2 -- balanced map"
      "\n       3 -- subtree-subset map"
      "\n       4 -- domain decomposition map"
      "\n    cutoff    -- cutoff for domain size w.r.t. # vertices"
      "\n       used only for the domain decomposition map"
      "\n    etreeFile -- output file, must be *.etreef or *.etreeb"
      "\n    mapFile   -- output file, must be *.ivf or *.ivb"
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
n1 = atoi(argv[3]) ;
n2 = atoi(argv[4]) ;
n3 = atoi(argv[5]) ;
maxzeros = atoi(argv[6]) ;
maxsize  = atoi(argv[7]) ;
nthread  = atoi(argv[8]) ;
maptype  = atof(argv[9]) ;
cutoff   = atof(argv[10]) ;
outETreeFileName = argv[11] ;
outMapIVfileName = argv[12] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n n1        -- %d" 
        "\n n2        -- %d" 
        "\n n3        -- %d" 
        "\n maxzeros  -- %d" 
        "\n maxsize   -- %d" 
        "\n nthread   -- %d" 
        "\n maptype   -- %d" 
        "\n cutoff    -- %f" 
        "\n etreeFile -- %s" 
        "\n mapFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, 
        maxzeros, maxsize, nthread, maptype, cutoff,
        outETreeFileName, outMapIVfileName) ;
fflush(msgFile) ;
if ( maptype < 0 || maptype > 4 ) {
   fprintf(stderr, "\n fatal error, maptype = %d, use "
           "\n 1 -- wrap map"
           "\n 2 -- balanced map"
           "\n 3 -- subtree-subset map"
           "\n 4 -- domain decomposition map\n", maptype) ;
   exit(-1) ;
}
/*
   ----------------------------
   create the grid graph object
   ----------------------------
*/
if ( n3 == 1 ) {
   egraph = EGraph_make9P(n1, n2, 1) ;
} else {
   egraph = EGraph_make27P(n1, n2, n3, 1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d grid EGraph", n1, n2, n3) ;
   EGraph_writeForHumanEye(egraph, msgFile) ;
   fflush(msgFile) ;
}
graph = EGraph_mkAdjGraph(egraph) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d grid Graph", n1, n2, n3) ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   get the nested dissection ordering
   ----------------------------------
*/
nvtx = n1*n2*n3 ;
newToOld = IVinit(nvtx, -1) ;
oldToNew = IVinit(nvtx, -1) ;
mkNDperm(n1, n2, n3, newToOld, 0, n1-1, 0, n2-1, 0, n3-1) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   oldToNew[newToOld[v]] = v ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n %d x %d x %d nd ordering", n1, n2, n3) ;
   IVfprintf(msgFile, nvtx, oldToNew) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------------
   create the fundamental supernode ETree object
   ---------------------------------------------
*/
nzerosIV = IV_new() ; 
IV_init(nzerosIV, nvtx, NULL) ;
IV_fill(nzerosIV, 0) ;
etree0 = ETree_new() ;
ETree_initFromGraphWithPerms(etree0, graph, newToOld, oldToNew) ;
etree1 = ETree_mergeFrontsOne(etree0, 0, nzerosIV) ;
nfronts[0] = ETree_nfront(etree1) ;
nfind[0]   = ETree_nFactorIndices(etree1) ;
nzf[0]     = ETree_nFactorEntries(etree1, SPOOLES_SYMMETRIC) ;
ops[0]     = ETree_nFactorOps(etree1, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n\n fs tree  : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[0], nfind[0], nzf[0], ops[0]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n fundamental elimination tree") ;
   ETree_writeForHumanEye(etree1, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   first step: try to absorb an only child
   ---------------------------------------
*/
etree2 = ETree_mergeFrontsOne(etree1, maxzeros, nzerosIV) ;
nfronts[1] = ETree_nfront(etree2) ;
nfind[1]   = ETree_nFactorIndices(etree2) ;
nzf[1]     = ETree_nFactorEntries(etree2, SPOOLES_SYMMETRIC) ;
ops[1]     = ETree_nFactorOps(etree2, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge 1  : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[1], nfind[1], nzf[1], ops[1]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after first merge") ;
   ETree_writeForHumanEye(etree2, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------------------
   second step: try to absorb all children
   ---------------------------------------
*/
etree3 = ETree_mergeFrontsAll(etree2, maxzeros, nzerosIV) ;
nfronts[2] = ETree_nfront(etree3) ;
nfind[2]   = ETree_nFactorIndices(etree3) ;
nzf[2]     = ETree_nFactorEntries(etree3, SPOOLES_SYMMETRIC) ;
ops[2]     = ETree_nFactorOps(etree3, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n merge 2  : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[2], nfind[2], nzf[2], ops[2]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after second merge") ;
   ETree_writeForHumanEye(etree3, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------
   third step: split the front tree
   --------------------------------
*/
etree4 = ETree_splitFronts(etree3, NULL, maxsize, 0) ;
nfronts[3] = ETree_nfront(etree4) ;
nfind[3]   = ETree_nFactorIndices(etree4) ;
nzf[3]     = ETree_nFactorEntries(etree4, SPOOLES_SYMMETRIC) ;
ops[3]     = ETree_nFactorOps(etree4, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n split    : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[3], nfind[3], nzf[3], ops[3]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree after split") ;
   ETree_writeForHumanEye(etree4, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   create the symbolic factorization object
   ----------------------------------------
*/
symbfacIVL = SymbFac_initFromGraph(etree4, graph) ;
nfronts[4] = ETree_nfront(etree4) ;
nfind[4]   = ETree_nFactorIndices(etree4) ;
nzf[4]     = ETree_nFactorEntries(etree4, SPOOLES_SYMMETRIC) ;
ops[4]     = ETree_nFactorOps(etree4, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
fprintf(msgFile,
        "\n final    : %8d fronts, %8d indices, %8d |L|, %12.0f ops",
        nfronts[4], nfind[4], nzf[4], ops[4]) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after symbolic factorization") ;
   ETree_writeForHumanEye(etree4, msgFile) ;
   fprintf(msgFile, "\n\n after symbolic factorization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
}
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree4, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from ETree_writeToFile(%p,%s)",
              rc, etree4, outETreeFileName) ;
   }
}
/*
   ----------------------------
   get the owners map IV object
   ----------------------------
*/
cumopsDV = DV_new() ;
DV_init(cumopsDV, nthread, NULL) ;
DV_fill(cumopsDV, 0.0) ;
switch ( maptype ) {
case 1 : /* wrap map */
   ownersIV = ETree_wrapMap(etree4, SPOOLES_REAL,
                            SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 2 : /* balanced map */
   ownersIV = ETree_balancedMap(etree4, SPOOLES_REAL,
                                SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 3 : /* subtree-subset map */
   ownersIV = ETree_subtreeSubsetMap(etree4, SPOOLES_REAL,
                                     SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 4 : /* dd map */
   msIV = ETree_msByNvtxCutoff(etree4, cutoff) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n multisector IV object") ;
      IV_writeForHumanEye(msIV, msgFile) ;
      fflush(msgFile) ;
   }
   ownersIV = ETree_ddMapNew(etree4, SPOOLES_REAL,
                             SPOOLES_SYMMETRIC, msIV, cumopsDV) ;
   IV_free(msIV) ;
   break ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n totalOps = %.0f", DV_sum(cumopsDV)) ;
   DVscale(DV_size(cumopsDV), DV_entries(cumopsDV),
            nthread/DV_sum(cumopsDV)) ;
   fprintf(msgFile, "\n\n cumopsDV") ;
   DV_writeForHumanEye(cumopsDV, msgFile) ;
   fflush(msgFile) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n ownersIV") ;
   IV_writeForHumanEye(ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the map IV object
   ---------------------------
*/
if ( strcmp(outMapIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(ownersIV, outMapIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write owners map to file %s",
           t2 - t1, outMapIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from IV_writeToFile(%p,%s)",
              rc, ownersIV, outMapIVfileName) ;
   }
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree0) ;
ETree_free(etree1) ;
ETree_free(etree2) ;
ETree_free(etree3) ;
ETree_free(etree4) ;
EGraph_free(egraph) ;
Graph_free(graph) ;
IVfree(newToOld) ;
IVfree(oldToNew) ;
IV_free(nzerosIV) ;
IV_free(ownersIV) ;
IVL_free(symbfacIVL) ;
DV_free(cumopsDV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
