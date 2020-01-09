/*  mkGridEgraph.c  */

#include "../EGraph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   generate an Egraph object for a 9pt or 27pt grid
   optionally write out to a file

   created -- 95nov08, cca
   -------------------------------------------------
*/
{
double   t1, t2 ;
int      msglvl, n1, n2, n3, ncomp, rc ;
EGraph   *eg ;
FILE     *msgFile ;

if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile n1 n2 n3 ncomp outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    n1       -- # of points in first grid direction"
      "\n    n2       -- # of points in second grid direction"
      "\n    n3       -- # of points in third grid direction"
      "\n    ncomp    -- # of components per grid point"
      "\n    outFile  -- output file, must be *.egraphf or *.egraphb"
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
n1    = atoi(argv[3]) ;
n2    = atoi(argv[4]) ;
n3    = atoi(argv[5]) ;
ncomp = atoi(argv[6]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n n1       -- %d" 
        "\n n2       -- %d" 
        "\n n3       -- %d" 
        "\n ncomp    -- %d" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], n1, n2, n3, ncomp, argv[7]) ;
fflush(msgFile) ;
/*
   ---------------------------
   create in the EGraph object
   ---------------------------
*/
MARKTIME(t1) ;
if ( n3 == 1 ) {
   eg = EGraph_make9P(n1, n2, ncomp) ;
} else {
   eg = EGraph_make27P(n1, n2, n3, ncomp) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : create the Egraph object", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n EGraph object") ;
   EGraph_writeForHumanEye(eg, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the EGraph object
   ---------------------------
*/
if ( strcmp(argv[7], "none") != 0 ) {
   MARKTIME(t1) ;
   rc = EGraph_writeToFile(eg, argv[7]) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write graph to file %s",
           t2 - t1, argv[7]) ;
   if ( rc != 1 ) {
      fprintf(msgFile, 
              "\n return value %d from Graph_writeToFile(%p,%s)",
              rc, eg, argv[7]) ;
   }
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
