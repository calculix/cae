/*  testMaps.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------
   test the map routines

   created -- 97jan15, cca
   -----------------------
*/
{
char     *inETreeFileName, *outIVfileName ;
double   cutoff, t1, t2 ;
DV       *cumopsDV ;
int      msglvl, nthread, rc, type ;
IV       *ownersIV ;
ETree    *etree ;
FILE     *msgFile ;

if ( argc != 8 ) {
   fprintf(stdout, 
    "\n\n usage : %s msglvl msgFile inETreeFile outIVfile "
    "\n         nthread type cutoff"
    "\n    msglvl      -- message level"
    "\n    msgFile     -- message file"
    "\n    inETreeFile -- input file, must be *.etreef or *.etreeb"
    "\n    outIVfile   -- output file, must be *.ivf or *.ivb"
    "\n    nthread     -- number of threads"
    "\n    type        -- type of map"
    "\n       1 -- wrap map"
    "\n       2 -- balanced map via a post-order traversal"
    "\n       3 -- subtree-subset map"
    "\n       4 -- old domain decomposition map"
    "\n       5 -- new domain decomposition map"
    "\n    cutoff -- cutoff used for the domain decomposition map"
    "\n       0 <= cutoff <= 1 used to define the multisector"
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
outIVfileName = argv[4] ;
nthread = atoi(argv[5]) ;
type    = atoi(argv[6]) ;
cutoff  = atof(argv[7]) ;
if ( type < 1 || type > 5 ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n type = %d, must be 1, 2 or 3",
           argv[0], type) ;
   exit(-1) ;
}
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inETreeFile -- %s" 
        "\n outIVfile   -- %s" 
        "\n nthread     -- %d" 
        "\n type        -- %d" 
        "\n cutoff      -- %f" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, outIVfileName,
        nthread, type, cutoff) ;
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
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------------------------------
   initialize the cumulative operations metric object
   --------------------------------------------------
*/
cumopsDV = DV_new() ;
DV_init(cumopsDV, nthread, NULL) ;
DV_fill(cumopsDV, 0.0) ;
/*
   -------------------------------
   create the owners map IV object
   -------------------------------
*/
switch ( type ) {
case 1 :
   ownersIV = ETree_wrapMap(etree, SPOOLES_REAL, 
                            SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 2 :
   ownersIV = ETree_balancedMap(etree, SPOOLES_REAL,
                                SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 3 :
   ownersIV = ETree_subtreeSubsetMap(etree, SPOOLES_REAL,
                                     SPOOLES_SYMMETRIC, cumopsDV) ;
   break ;
case 4 :
   ownersIV = ETree_ddMap(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC, 
                          cumopsDV, cutoff) ;
   break ;
case 5 :
   {  IV  *msIV ;
/*
   msIV = ETree_msByNopsCutoff(etree, cutoff, 1) ;
*/
   msIV = ETree_msByNvtxCutoff(etree, cutoff) ;
   ownersIV = ETree_ddMapNew(etree, SPOOLES_REAL,
                             SPOOLES_SYMMETRIC, msIV, cumopsDV) ;
   IV_free(msIV) ;
   }
   break ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n totalOps = %.0f", DV_sum(cumopsDV)) ;
   DVscale(DV_size(cumopsDV), DV_entries(cumopsDV),
            nthread/DV_sum(cumopsDV)) ;
   fprintf(msgFile, "\n\n cumopsDV") ;
   DV_writeForHumanEye(cumopsDV, msgFile) ;
   fprintf(msgFile, "\n\n ownersIV") ;
   IV_writeForHumanEye(ownersIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------
   write out the IV object
   --------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(ownersIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write ownersIV to file %s",
           t2 - t1, outIVfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
           rc, ownersIV, outIVfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
DV_free(cumopsDV) ;
IV_free(ownersIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
