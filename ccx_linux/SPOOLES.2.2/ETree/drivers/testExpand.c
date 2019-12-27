/*  testExpand.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   read in an ETree object and an equivalence map,
   expand the ETree object and optionally write to a file.

   created -- 98sep05, cca
   -------------------------------------------------------
*/
{
char     *inEqmapFileName, *inETreeFileName, *outETreeFileName ;
double   t1, t2 ;
ETree    *etree, *etree2 ;
FILE     *msgFile ;
int      msglvl, rc ;
IV       *eqmapIV ;

if ( argc != 6 ) {
   fprintf(stdout, 
   "\n\n usage : %s msglvl msgFile inETreeFile inEqmapFile outETreeFile"
   "\n    msglvl       -- message level"
   "\n    msgFile      -- message file"
   "\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
   "\n    inEqmapFile  -- input file, must be *.ivf or *.ivb"
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
inETreeFileName  = argv[3] ;
inEqmapFileName  = argv[4] ;
outETreeFileName = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl       -- %d" 
        "\n msgFile      -- %s" 
        "\n inETreeFile  -- %s" 
        "\n inEqmapFile  -- %s" 
        "\n outETreeFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inETreeFileName, inEqmapFileName, outETreeFileName) ;
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
   -----------------------
   expand the ETree object
   -----------------------
*/
etree2 = ETree_expand(etree, eqmapIV) ;
fprintf(msgFile, "\n\n after expanding the ETree object") ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree2, msgFile) ;
} else {
   ETree_writeStats(etree2, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree2, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, etree2, outETreeFileName) ;
}
/*
   ---------------------
   free the ETree object
   ---------------------
*/
ETree_free(etree) ;
IV_free(eqmapIV) ;
ETree_free(etree2) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
