/*  testDM.c  */

#include "../BPG.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   read BPG from file and get the Dulmage-Mendelsohn decomposition

   created -- 96mar08, cca
   ---------------------------------------------------------------
*/
{
char     *inBPGFileName ;
double   t1, t2 ;
int      ierr, msglvl, rc ;
int      *dmflags, *stats ;
BPG      *bpg ;
FILE     *msgFile ;

if ( argc != 4 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile "
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.bpgf or *.bpgb"
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
inBPGFileName  = argv[3] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inBPGFileName) ;
fflush(msgFile) ;
/*
   ----------------------
   read in the BPG object
   ----------------------
*/
if ( strcmp(inBPGFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
bpg = BPG_new() ;
MARKTIME(t1) ;
rc = BPG_readFromFile(bpg, inBPGFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inBPGFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from BPG_readFromFile(%p,%s)",
           rc, bpg, inBPGFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading BPG object from file %s",
        inBPGFileName) ;
if ( msglvl > 2 ) {
   BPG_writeForHumanEye(bpg, msgFile) ;
} else {
   BPG_writeStats(bpg, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------
   test out the max flow DMdecomposition method
   --------------------------------------------
*/
dmflags = IVinit(bpg->nX + bpg->nY, -1) ;
stats   = IVinit(6, 0) ;
MARKTIME(t1) ;
BPG_DMviaMaxFlow(bpg, dmflags, stats, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %9.5f : find DM via maxflow", t2 - t1) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n\n BPG_DMviaMaxFlow"
           "\n |X_I| = %6d, |X_E| = %6d, |X_R| = %6d"
           "\n |Y_I| = %6d, |Y_E| = %6d, |Y_R| = %6d",
           stats[0], stats[1], stats[2],
           stats[3], stats[4], stats[5]) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n dmflags") ;
   IVfp80(msgFile, bpg->nX + bpg->nY, dmflags, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ------------------------------------------
   test out the matching DMcomposition method
   ------------------------------------------
*/
IVfill(bpg->nX + bpg->nY, dmflags, -1) ;
IVfill(6, stats, -1) ;
MARKTIME(t1) ;
BPG_DMdecomposition(bpg, dmflags, stats, msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n\n CPU %9.5f : find DM via matching", t2 - t1) ;
if ( msglvl > 0 ) {
   fprintf(msgFile, 
           "\n\n BPG_DMdecomposition"
           "\n |X_I| = %6d, |X_E| = %6d, |X_R| = %6d"
           "\n |Y_I| = %6d, |Y_E| = %6d, |Y_R| = %6d",
           stats[0], stats[1], stats[2],
           stats[3], stats[4], stats[5]) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n dmflags") ;
   IVfp80(msgFile, bpg->nX + bpg->nY, dmflags, 80, &ierr) ;
   fflush(msgFile) ;
}
/*
   ----------------
   free the storage
   ----------------
*/
IVfree(dmflags) ;
IVfree(stats) ;
BPG_free(bpg) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
