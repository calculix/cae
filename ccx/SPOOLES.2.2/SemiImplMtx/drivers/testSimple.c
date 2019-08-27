/*  testSimple.c  */

#include "../SemiImplMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   (1) read in a FrontMtx from a file,
   (2) read in an IV object that has the domain/schur complement
       information for the fronts
   (3) create the SemiImplMtx object

   created -- 98oct16, cca
   -------------------------------------------------
*/
{
char          *inFrontMtxFileName, *inInpMtxFileName, *inIVfileName ;
double        t1, t2 ;
FILE          *msgFile ;
FrontMtx      *frontmtx ;
InpMtx        *inpmtx ;
int           msglvl, rc ;
IV            *frontmapIV ;
SemiImplMtx   *semimtx ;

if ( argc != 6 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inFrontMtxFile inInpMtxFile inIVfile"
"\n    msglvl         -- message level"
"\n    msgFile        -- message file"
"\n    inFrontMtxFile -- input file, must be *.frontmtxf or *.frontmtxb"
"\n    inInpMtxFile   -- input file, must be *.inpmtxf or *.inpmtxb"
"\n    inIVfile       -- input file, must be *.ivf or *.ivb"
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
inFrontMtxFileName = argv[3] ;
inInpMtxFileName   = argv[4] ;
inIVfileName       = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl         -- %d" 
        "\n msgFile        -- %s" 
        "\n inFrontMtxFile -- %s" 
        "\n inInpMtxFile   -- %s" 
        "\n inIVfile       -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inFrontMtxFileName, 
        inInpMtxFileName, inIVfileName) ;
fflush(msgFile) ;
/*
   ---------------------------
   read in the FrontMtx object
   ---------------------------
*/
if ( strcmp(inFrontMtxFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
frontmtx = FrontMtx_new() ;
MARKTIME(t1) ;
rc = FrontMtx_readFromFile(frontmtx, inFrontMtxFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in frontmtx from file %s",
        t2 - t1, inFrontMtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from FrontMtx_readFromFile(%p,%s)",
           rc, frontmtx, inFrontMtxFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading FrontMtx object from file %s",
        inFrontMtxFileName) ;
if ( msglvl > 2 ) {
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
} else {
   FrontMtx_writeStats(frontmtx, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------
   read in the InpMtx object
   ---------------------------
*/
if ( strcmp(inInpMtxFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
inpmtx = InpMtx_new() ;
MARKTIME(t1) ;
rc = InpMtx_readFromFile(inpmtx, inInpMtxFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in inpmtx from file %s",
        t2 - t1, inInpMtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from InpMtx_readFromFile(%p,%s)",
           rc, inpmtx, inInpMtxFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
        inInpMtxFileName) ;
if ( msglvl > 2 ) {
   InpMtx_writeForHumanEye(inpmtx, msgFile) ;
} else {
   InpMtx_writeStats(inpmtx, msgFile) ;
}
fflush(msgFile) ;
/*
   ---------------------------------
   read in the fronts' map IV object
   ---------------------------------
*/
if ( strcmp(inIVfileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
frontmapIV = IV_new() ;
MARKTIME(t1) ;
rc = IV_readFromFile(frontmapIV, inIVfileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in frontmapIV from file %s",
        t2 - t1, inIVfileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from IV_readFromFile(%p,%s)",
           rc, frontmapIV, inIVfileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading IV object from file %s",
        inIVfileName) ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(frontmapIV, msgFile) ;
} else {
   IV_writeStats(frontmapIV, msgFile) ;
}
fflush(msgFile) ;
/*
   -------------------------------
   create the semi-implicit matrix
   -------------------------------
*/
semimtx = SemiImplMtx_new() ;
rc = SemiImplMtx_initFromFrontMtx(semimtx, frontmtx, inpmtx, 
                                  frontmapIV, msglvl, msgFile) ;
fprintf(msgFile, "\n\n Semi-implicit matrix") ;
SemiImplMtx_writeForHumanEye(semimtx, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
{
ETree   *etree ;
IVL     *symbfacIVL ;

etree = frontmtx->frontETree ;
symbfacIVL = frontmtx->symbfacIVL ;
FrontMtx_free(frontmtx) ;
ETree_free(etree) ;
IVL_free(symbfacIVL) ;
}
SemiImplMtx_free(semimtx) ;
InpMtx_free(inpmtx) ;
IV_free(frontmapIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
