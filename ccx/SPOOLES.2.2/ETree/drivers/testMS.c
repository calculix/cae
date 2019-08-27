/*  testMS.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------
   generate multisectors

   created -- 96oct31, cca
   -----------------------
*/
{
char     *inETreeFileName, *outIVfileName ;
double   cutoff, domops, totops, t1, t2 ;
double   *opsreg ;
DV       *opsDV ;
int      domnzf, domnvtx, flag, depth, ireg, msglvl, nreg, rc,
         totnvtx, totnzf ;
int      *nvtxreg, *nzfreg ;
IV       *msIV, *nvtxIV, *nzfIV ;
ETree    *etree ;
FILE     *msgFile ;

if ( argc != 7 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile outIVfile flag cutoff"
"\n    msglvl      -- message level"
"\n    msgFile     -- message file"
"\n    inETreeFile -- input file, must be *.etreef or *.etreeb"
"\n    outIVfile   -- output file, must be *.ivf or *.ivb"
"\n    flag        -- type of multisector"
"\n       flag = 1 --> multisector via depth of front"
"\n       flag = 2 --> multisector via # of vertices in subtree"
"\n       flag = 3 --> multisector via # of entries in subtree"
"\n       flag = 4 --> multisector via # of operations in subtree"
"\n    cutoff      -- cutoff point for multisector"
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
inETreeFileName = argv[3] ;
outIVfileName   = argv[4] ;
flag            = atoi(argv[5]) ;
if ( flag == 1 ) {
   depth = atoi(argv[6]) ;
   fprintf(msgFile, 
           "\n %s "
           "\n msglvl      -- %d" 
           "\n msgFile     -- %s" 
           "\n inETreeFile -- %s" 
           "\n outIVfile   -- %s" 
           "\n flag        -- %d" 
           "\n depth       -- %d" 
           "\n",
           argv[0], msglvl, argv[2], inETreeFileName, outIVfileName,
           flag, depth) ;
   fflush(msgFile) ;
} else if ( 2 <= flag && flag <= 4 ) {
   cutoff = atof(argv[6]) ;
   fprintf(msgFile, 
           "\n %s "
           "\n msglvl      -- %d" 
           "\n msgFile     -- %s" 
           "\n inETreeFile -- %s" 
           "\n outIVfile   -- %s" 
           "\n flag        -- %d" 
           "\n cutoff      -- %f" 
           "\n",
           argv[0], msglvl, argv[2], inETreeFileName, outIVfileName,
           flag, cutoff) ;
   fflush(msgFile) ;
}
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
   --------------------
   make the multisector
   --------------------
*/
switch ( flag ) {
case  1 :
   msIV = ETree_msByDepth(etree, depth) ;
   break ;
case  2 :
   msIV = ETree_msByNvtxCutoff(etree, cutoff) ;
   break ;
case  3 :
   msIV = ETree_msByNentCutoff(etree, cutoff, SPOOLES_SYMMETRIC) ;
   break ;
case  4 :
   msIV = ETree_msByNopsCutoff(etree, cutoff, 
                               SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   break ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n msIV") ;
   IV_writeForHumanEye(msIV, msgFile) ;
}
/*
   -----------------------
   generate the statistics
   -----------------------
*/
nvtxIV = IV_new() ;
nzfIV  = IV_new() ;
opsDV  = DV_new() ;
ETree_msStats(etree, msIV, nvtxIV, nzfIV, opsDV, 
              SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
nreg = IV_size(nvtxIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n msIV") ;
   IV_writeForHumanEye(msIV, msgFile) ;
}
nvtxreg = IV_entries(nvtxIV) ;
nzfreg  = IV_entries(nzfIV) ;
opsreg  = DV_entries(opsDV) ;
domnvtx = IVsum(nreg-1, nvtxreg+1) ;
domnzf  = IVsum(nreg-1, nzfreg+1) ;
domops  = DVsum(nreg-1, opsreg+1) ;
fprintf(msgFile, 
"\n region  vertices    entries   operations   metric/(avg domain)") ;
for ( ireg = 0 ; ireg < nreg ; ireg++ ) {
   fprintf(msgFile, "\n %5d %10d %10d %12.0f %6.3f %6.3f %6.3f",
           ireg, nvtxreg[ireg], nzfreg[ireg], opsreg[ireg],
           ((double) nvtxreg[ireg]*(nreg-1))/domnvtx,
           ((double) nzfreg[ireg]*(nreg-1))/domnzf,
           opsreg[ireg]*(nreg-1)/domops) ;
}
totnvtx = IV_sum(nvtxIV) ;
totnzf  = IV_sum(nzfIV) ;
totops  = DV_sum(opsDV) ;
fprintf(msgFile, 
"\n\n                    nvtx   %%         nzf   %%           ops    %%"
        "\n domains          %6d %5.2f %9d %5.2f %12.0f %5.3f"
        "\n schur complement %6d %5.2f %9d %5.2f %12.0f %5.3f"
        "\n total            %6d       %9d       %12.0f      ",
        domnvtx, (100.*domnvtx)/totnvtx,
        domnzf, (100.*domnzf)/totnzf,
        domops, (100.*domops)/totops,
        totnvtx - domnvtx, (100.*(totnvtx - domnvtx))/totnvtx,
        totnzf - domnzf, (100.*(totnzf - domnzf))/totnzf,
        totops - domops, (100.*(totops - domops))/totops,
        totnvtx, totnzf, totops) ;

fprintf(msgFile, "\n\n FrontETree : %d entries, %.0f ops",
        ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC),
        ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC)) ;
/*
   -----------------------
   write out the IV object
   -----------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(msIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write msIV to file %s",
           t2 - t1, outIVfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
           rc, msIV, outIVfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
IV_free(msIV) ;
IV_free(nvtxIV) ;
IV_free(nzfIV) ;
DV_free(opsDV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
