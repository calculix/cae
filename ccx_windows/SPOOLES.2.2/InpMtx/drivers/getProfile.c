/*  getProfile.c  */

#include "../InpMtx.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------
   get a log10 profile of the matrix entries

   created -- 97feb14, cca
   -----------------------------------------
*/
{
int       ii, msglvl, nbig, npts, nsmall, nzero, rc, size ;
InpMtx   *inpmtx ;
double    taubig, tausmall ;
double    *xvec, *yvec ;
DV        *xDV, *yDV ;
FILE      *msgFile ;

if ( argc != 7 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile npts tausmall taubig"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.dinpmtxf or *.dinpmtxb"
      "\n    npts     -- number of points in the profile curve"
      "\n    tausmall -- lower cutoff"
      "\n    taubig   -- upper cutoff"
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
npts     = atoi(argv[4]) ;
tausmall = atof(argv[5]) ;
taubig   = atof(argv[6]) ;
fprintf(msgFile, 
        "\n %% %s "
        "\n %% msglvl   -- %d" 
        "\n %% msgFile  -- %s" 
        "\n %% inFile   -- %s" 
        "\n %% npts     -- %d" 
        "\n %% tausmall -- %e" 
        "\n %% taubig   -- %e" 
        "\n",
        argv[0], msglvl, argv[2], argv[3], npts, tausmall, taubig) ;
fflush(msgFile) ;
/*
   --------------------------
   read in the InpMtx object
   --------------------------
*/
if ( strcmp(argv[3], "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
inpmtx = InpMtx_new() ;
rc = InpMtx_readFromFile(inpmtx, argv[3]) ;
fprintf(msgFile, 
        "\n %% return value %d from InpMtx_readFromFile(%p,%s)",
        rc, inpmtx, argv[3]) ;
if ( rc != 1 ) {
   exit(-1) ;
}
fprintf(msgFile, "\n\n %d entries", inpmtx->nent) ;
/*
   ---------------
   get the profile
   ---------------
*/
xDV = DV_new() ;
yDV = DV_new() ;
InpMtx_log10profile(inpmtx, npts, xDV, yDV, tausmall, taubig,
                    &nzero, &nsmall, &nbig) ;
fprintf(msgFile, 
        "\n %% %8d zero entries "
        "\n %% %8d entries smaller than %20.12e in magnitude"
        "\n %% %8d entries larger  than %20.12e in magnitude",
        nzero, nsmall, tausmall, nbig, taubig) ;
DV_sizeAndEntries(xDV, &size, &xvec) ;
DV_sizeAndEntries(yDV, &size, &yvec) ;
fprintf(msgFile, "\n data = [ ...") ;
for ( ii = 0 ; ii < size ; ii++ ) {
   fprintf(msgFile, "\n %20.12e %20.12e", xvec[ii], yvec[ii]) ;
}
fprintf(msgFile, " ] ; ") ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
DV_free(xDV) ;
DV_free(yDV) ;
InpMtx_free(inpmtx) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
