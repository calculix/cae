/*  testIO.c  */

#include "../EGraph.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------
   test EGraph_readFromFile and EGraph_writeToFile,
   useful for translating between formatted *.egraphf
   and binary *.egraphb files.

   created -- 95nov03, cca
   --------------------------------------------------
*/
{
double    t1, t2 ;
int       msglvl, rc ;
EGraph    egraph ;
FILE      *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.egraphf or *.egraphb"
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
fprintf(msgFile, 
        "\n %s "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        argv[0], msglvl, argv[2], argv[3], argv[4]) ;
fflush(msgFile) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
EGraph_setDefaultFields(&egraph) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after setting default fields") ;
   EGraph_writeForHumanEye(&egraph, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the EGraph object
   -------------------------
*/
if ( strcmp(argv[3], "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = EGraph_readFromFile(&egraph, argv[3]) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in egraph from file %s",
        t2 - t1, argv[3]) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from EGraph_readFromFile(%p,%s)",
           rc, &egraph, argv[3]) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading EGraph object from file %s",
           argv[3]) ;
   EGraph_writeForHumanEye(&egraph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the EGraph object
   ---------------------------
*/
if ( strcmp(argv[4], "none") != 0 ) {
   MARKTIME(t1) ;
   rc = EGraph_writeToFile(&egraph, argv[4]) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write egraph to file %s",
           t2 - t1, argv[4]) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from EGraph_writeToFile(%p,%s)",
           rc, &egraph, argv[4]) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
