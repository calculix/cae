/*  testExpand.c  */

#include "../IVL.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   read in an IVL object and an equivalence map,
   expand the IVL object and optionally write to a file.

   created -- 98sep05, cca
   -------------------------------------------------------
*/
{
char     *inEqmapFileName, *inIVLfileName, *outIVLfileName ;
double   t1, t2 ;
IVL      *ivl, *ivl2 ;
FILE     *msgFile ;
int      msglvl, rc ;
IV       *eqmapIV ;

if ( argc != 6 ) {
   fprintf(stdout, 
   "\n\n usage : %s msglvl msgFile inIVLfile inEqmapFile outIVLfile"
   "\n    msglvl      -- message level"
   "\n    msgFile     -- message file"
   "\n    inIVLfile   -- input file, must be *.ivlf or *.ivlb"
   "\n    inEqmapFile -- input file, must be *.ivf or *.ivb"
   "\n    outIVLfile  -- output file, must be *.ivlf or *.ivlb"
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
inIVLfileName   = argv[3] ;
inEqmapFileName = argv[4] ;
outIVLfileName  = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inIVLfile   -- %s" 
        "\n inEqmapFile -- %s" 
        "\n outIVLfile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inIVLfileName, inEqmapFileName, outIVLfileName) ;
fflush(msgFile) ;
/*
   ----------------------
   read in the IVL object
   ----------------------
*/
if ( strcmp(inIVLfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
ivl = IVL_new() ;
IVL_init1(ivl, IVL_CHUNKED, 0) ;
MARKTIME(t1) ;
rc = IVL_readFromFile(ivl, inIVLfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in ivl from file %s",
        t2 - t1, inIVLfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IVL_readFromFile(%p,%s)",
           rc, ivl, inIVLfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IVL object from file %s",
        inIVLfileName) ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(ivl, msgFile) ;
} else {
   IVL_writeStats(ivl, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------------------
   read in the equivalence map IV object
   -------------------------------------
*/
if ( strcmp(inEqmapFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
eqmapIV = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(eqmapIV, inEqmapFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in eqmapIV from file %s",
        t2 - t1, inEqmapFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
           rc, eqmapIV, inEqmapFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inEqmapFileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(eqmapIV, msgFile) ;
} else {
   IV_writeStats(eqmapIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------
   expand the IVL object
   ---------------------
*/
MARKTIME(t1) ;
ivl2 = IVL_expand(ivl, eqmapIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : expand the IVL object", t2 - t1) ;
fprintf(msgFile, "\n\n after expanding the IVL object") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(ivl2, msgFile) ;
} else {
   IVL_writeStats(ivl2, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------
   write out the IVL object
   ------------------------
*/
if ( strcmp(outIVLfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IVL_writeToFile(ivl2, outIVLfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write ivl to file %s",
           t2 - t1, outIVLfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IVL_writeToFile(%p,%s)",
           rc, ivl2, outIVLfileName) ;
}
/*
   ---------------------
   free the IVL object
   ---------------------
*/
IVL_free(ivl) ;
IV_free(eqmapIV) ;
IVL_free(ivl2) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
