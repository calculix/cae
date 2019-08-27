/*  extractBPG.c  */

#include "../BPG.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------------------
   extract a bipartite graph from a graph and a two-set partition

   created -- 96nov02, cca
   --------------------------------------------------------------
*/
{
char     *inGraphFileName, *inCompidsIVfileName, 
         *outBPGfileName, *outMapIVfileName ;
double   t1, t2 ;
int      msglvl, rc ;
BPG      *bpg ;
FILE     *msgFile ;
Graph    *graph ;
int      icomp, ierr, nvtx ;
int      *cmap, *compids, *indX, *indY ;
IV       *compidsIV ;

if ( argc != 8 ) {
   fprintf(stdout, 
   "\n\n usage : %s msglvl msgFile inGraphFile inCompidsIVfile"
   "\n         icomp outBPGfile"
   "\n    msglvl          -- message level"
   "\n    msgFile         -- message file"
   "\n    inGraphFile     -- input file, must be *.graphf or *.graphb"
   "\n    inCompidsIVfile -- input file, must be *.ivf or *.ivb"
   "\n    icomp           -- component for Y nodes"
   "\n    outMapIVfile    -- map output file, must be *.ivf or *.ivb"
   "\n    outBPGfile      -- output file, must be *.bpgf or *.bpgb"
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
inGraphFileName     = argv[3] ;
inCompidsIVfileName = argv[4] ;
icomp               = atoi(argv[5]) ;
outMapIVfileName    = argv[6] ;
outBPGfileName      = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl          -- %d" 
        "\n msgFile         -- %s" 
        "\n inGraphFile     -- %s" 
        "\n inCompidsIVfile -- %s" 
        "\n icomp           -- %d" 
        "\n outMapIVfile    -- %s" 
        "\n outBPGfile      -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inGraphFileName, inCompidsIVfileName,
        icomp, outMapIVfileName, outBPGfileName) ;
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
nvtx = graph->nvtx ;
/*
   ------------------------
   read in the IV object
   ------------------------
*/
if ( strcmp(inCompidsIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
compidsIV = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(compidsIV, inCompidsIVfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in compidsIV from file %s",
        t2 - t1, inCompidsIVfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, compidsIV, inCompidsIVfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inCompidsIVfileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(compidsIV, msgFile) ;
} else {
   IV_writeStats(compidsIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------------------------------
   extract out the bipartite graph that corresponds to
   the separator and its boundary in component icomp
   ---------------------------------------------------
*/
compids = IV_entries(compidsIV) ;
cmap = IVinit(nvtx, -1) ;
indX = IVinit(nvtx, -1) ;
indY = IVinit(nvtx, -1) ;
bpg = BPG_new() ;
BPG_initFromColoring(bpg, graph, compids, 0, icomp, cmap, indX, indY) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n cmap[]") ;
   IVfp80(msgFile, nvtx, cmap, 80, &ierr) ;
   fprintf(msgFile, "\n\n indX[]") ;
   IVfp80(msgFile, bpg->nX, indX, 80, &ierr) ;
   fprintf(msgFile, "\n\n indY[]") ;
   IVfp80(msgFile, bpg->nY, indY, 80, &ierr) ;
   fprintf(msgFile, "\n\n bipartite graph") ;
   BPG_writeForHumanEye(bpg, msgFile) ;
} else {
   BPG_writeStats(bpg, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------------------------------------------------
   if the map vector is requested, create it and write it to a file
   ----------------------------------------------------------------
*/
if ( strcmp(outMapIVfileName, "none") != 0 ) {
   IV   mapIV ;

   IV_setDefaultFields(&mapIV) ;
   IV_init(&mapIV, bpg->nX + bpg->nY, NULL) ;
   IVcopy(bpg->nX, IV_entries(&mapIV), indX) ;
   IVcopy(bpg->nY, IV_entries(&mapIV) + bpg->nX, indY) ;
   if ( msglvl > 0 ) {
      IV_writeForHumanEye(&mapIV, msgFile) ;
      fflush(msgFile) ;
   }
   MARKTIME(t1) ;
   rc = IV_writeToFile(&mapIV, outMapIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write mapIV to file %s",
           t2 - t1, outMapIVfileName) ;
}
/*
   ------------------------
   write out the BPG object
   ------------------------
*/
if ( strcmp(outBPGfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = BPG_writeToFile(bpg, outBPGfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, outBPGfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from BPG_writeToFile(%p,%s)",
           rc, bpg, outBPGfileName) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
