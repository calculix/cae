/*  testMemory.c  */

#include "../misc.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------
   test the memory subsystem

   created -- 98oct16, cca
   ----------------------------------------
*/
{
char     **p_chunks ;
double   t1, t2 ;
double   *cpus ;
int      chunksize, ichunk, jchunk, msglvl, nchunk ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
     "\n\n usage : %s msglvl msgFile nchunk chunksize"
     "\n    msglvl    -- message level"
     "\n    msgFile   -- message file"
     "\n    nchunk    -- # of chunks"
     "\n    chunksize -- size of each chunk in bytes"
     "\n",
argv[0]) ;
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
nchunk    = atoi(argv[3]) ;
chunksize = atoi(argv[4]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl    -- %d" 
        "\n msgFile   -- %s" 
        "\n nchunk    -- %d" 
        "\n chunksize -- %d" 
        "\n",
        argv[0], msglvl, argv[2], nchunk, chunksize) ;
fflush(msgFile) ;
/*
   ---------------------------------------
   allocate the pointers and cpus[] arrays
   ---------------------------------------
*/
p_chunks = (char **) malloc(nchunk*sizeof(char *)) ;
if ( p_chunks == NULL ) {
   fprintf(msgFile, "\n unable to allocate p_chunks") ;
   fflush(msgFile) ;
   return(-1) ;
}
cpus = (double *) malloc(2*nchunk*sizeof(double)) ;
if ( cpus == NULL ) {
   fprintf(msgFile, "\n unable to allocate cpus") ;
   fflush(msgFile) ;
   return(-1) ;
}
for ( ichunk = 0 ; ichunk < nchunk ; ichunk++ ) {
   p_chunks[ichunk] = NULL ;
   cpus[2*ichunk]   = cpus[2*ichunk+1] = 0.0 ;
}
/*
   -------------------
   allocate the chunks
   -------------------
*/
for ( ichunk = 0 ; ichunk < nchunk ; ichunk++ ) {
   MARKTIME(t1) ;
   p_chunks[ichunk] = (char *) malloc(chunksize * sizeof(char)) ;
   MARKTIME(t2) ;
   cpus[2*ichunk] = t2 - t1 ;
   if ( p_chunks[ichunk] == NULL ) {
      fprintf(msgFile, "\n unable to allocate p_chunks[%d]", ichunk) ;
      fflush(msgFile) ;
      break ;
   }
}
fprintf(msgFile, "\n\n %d chunks of size %d bytes allocated",
       ichunk, chunksize) ;
fprintf(msgFile, "\n %d total bytes allocated", ichunk*chunksize) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
for ( jchunk = 0 ; jchunk < ichunk ; jchunk++ ) {
   MARKTIME(t1) ;
   free(p_chunks[jchunk]) ;
   MARKTIME(t2) ;
   cpus[2*jchunk+1] = t2 - t1 ;
}
/*
   --------------------
   print the statistics
   --------------------
*/
fprintf(msgFile, "\n                    CPU TIME "
        "\n    chunk      malloc       free") ;
for ( jchunk = 0 ; jchunk < ichunk ; jchunk++ ) {
   fprintf(msgFile, "\n %8d    %8.4f    %8.4f",
           ichunk, cpus[2*ichunk], cpus[2*ichunk+1]) ;
}
free(p_chunks) ;
free(cpus) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
