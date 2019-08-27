/*  testIO.c  */

#include "../SubMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ------------------------------------------------------
   test SubMtx_readFromFile and SubMtx_writeToFile,
   useful for translating between formatted *.submtxf
   and binary *.submtxb files.

   created -- 98may04, cca
   ------------------------------------------------------
*/
{
char       *inSubMtxFileName, *matlabFileName, *outSubMtxFileName ;
double     t1, t2 ;
int        msglvl, rc ;
SubMtx     *mtx ;
FILE       *fp, *msgFile ;

if ( argc != 6 ) {
   fprintf(stdout, 
     "\n\n usage : %s msglvl msgFile inFile outFile matlabFile"
     "\n    msglvl     -- message level"
     "\n    msgFile    -- message file"
     "\n    inFile     -- input file, must be *.dmtxf or *.dmtxb"
     "\n    outFile    -- output file, must be *.dmtxf or *.dmtxb"
     "\n    matlabFile -- output file, must be *.m"
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
inSubMtxFileName  = argv[3] ;
outSubMtxFileName = argv[4] ;
matlabFileName      = argv[5] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl     -- %d" 
        "\n msgFile    -- %s" 
        "\n inFile     -- %s" 
        "\n outFile    -- %s" 
        "\n matlabFile -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inSubMtxFileName, outSubMtxFileName, matlabFileName) ;
fflush(msgFile) ;
/*
   ----------------------------
   read in the SubMtx object
   ----------------------------
*/
if ( strcmp(inSubMtxFileName, "none") == 0 ) { fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
mtx = SubMtx_new() ;
MARKTIME(t1) ;
rc = SubMtx_readFromFile(mtx, inSubMtxFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in dmtx from file %s",
        t2 - t1, inSubMtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from SubMtx_readFromFile(%p,%s)",
           rc, mtx, inSubMtxFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading SubMtx object from file %s",
        inSubMtxFileName) ;
fflush(msgFile) ;
if ( msglvl > 2 ) {
   SubMtx_writeForHumanEye(mtx, msgFile) ;
} else {
   SubMtx_writeStats(mtx, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------------------------------
   hack to convert from row major to column major
   ----------------------------------------------
*/
/*
{
SubMtx   *mtx2 ;
double   *ent1, *ent2, *pXij, *pYij ;
int      *colind, *colind2, *rowind, *rowind2 ;
int      inc1, inc2, irow, jcol, ncol, nrow ;

mtx2 = SubMtx_new() ;
SubMtx_init(mtx2, DMTX_DENSE_COLUMNS, 0, 0, 
          mtx->nrow, mtx->ncol, mtx->nent) ;
SubMtx_denseInfo(mtx, &nrow, &ncol, &inc1, &inc2, &ent1) ;
SubMtx_denseInfo(mtx2, &nrow, &ncol, &inc1, &inc2, &ent2) ;
SubMtx_rowIndices(mtx, &nrow, &rowind) ;
SubMtx_rowIndices(mtx2, &nrow, &rowind2) ;
IVcopy(nrow, rowind2, rowind) ;
SubMtx_columnIndices(mtx, &ncol, &colind) ;
SubMtx_columnIndices(mtx2, &ncol, &colind2) ;
IVcopy(ncol, colind2, colind) ;
for ( jcol = 0 ; jcol < ncol ; jcol++ ) {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      pXij = SubMtx_locationOfEntry(mtx,  irow, jcol) ;
      pYij = SubMtx_locationOfEntry(mtx2, irow, jcol) ;
      *pYij = *pXij ;
   }
}
SubMtx_free(mtx) ;
mtx = mtx2 ;
}
*/
/*
   ----------------------------
   write out the SubMtx object
   ----------------------------
*/
if ( strcmp(outSubMtxFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = SubMtx_writeToFile(mtx, outSubMtxFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write mtx to file %s",
           t2 - t1, outSubMtxFileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from SubMtx_writeToFile(%p,%s)",
           rc, mtx, outSubMtxFileName) ;
}
/*
   ----------------------------------------------
   write out the SubMtx object to a matlab file
   ----------------------------------------------
*/
if ( strcmp(matlabFileName, "none") != 0 ) {
   if ( (fp = fopen(matlabFileName, "a")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], matlabFileName) ;
      return(-1) ;
   }
   MARKTIME(t1) ;
   SubMtx_writeForMatlab(mtx, "rhs", fp) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write mtx to file %s",
           t2 - t1, matlabFileName) ;
}
/*
   -------------------------
   free the SubMtx object
   -------------------------
*/
SubMtx_free(mtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
