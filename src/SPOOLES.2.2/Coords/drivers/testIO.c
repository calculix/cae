/*  testIO.c  */

#include "../Coords.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   --------------------------------------------------
   test Coords_readFromFile and Coords_writeToFile,
   useful for translating between formatted *.coordsf
   and binary *.coordsb files.

   created -- 95dec17, cca
   --------------------------------------------------
*/
{
int      msglvl, rc ;
Coords   coords ;
FILE     *msgFile ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : testIO msglvl msgFile inFile outFile"
      "\n    msglvl   -- message level"
      "\n    msgFile  -- message file"
      "\n    inFile   -- input file, must be *.coordsf or *.coordsb"
      "\n    outFile  -- output file, must be *.coordsf or *.coordsb"
      "\n") ;
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
        "\n testIO "
        "\n msglvl   -- %d" 
        "\n msgFile  -- %s" 
        "\n inFile   -- %s" 
        "\n outFile  -- %s" 
        "\n",
        msglvl, argv[2], argv[3], argv[4]) ;
fflush(msgFile) ;
/*
   ----------------------
   set the default fields
   ----------------------
*/
Coords_setDefaultFields(&coords) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after setting default fields") ;
   Coords_writeForHumanEye(&coords, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the Coords object
   -------------------------
*/
if ( strcmp(argv[3], "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from\n") ;
   exit(0) ;
}
rc = Coords_readFromFile(&coords, argv[3]) ;
fprintf(msgFile, "\n return value %d from Coords_readFromFile(%p,%s)",
        rc, &coords, argv[3]) ;
if ( rc != 1 ) {
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading Coords object from file %s",
           argv[3]) ;
   Coords_writeForHumanEye(&coords, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the Coords object
   ---------------------------
*/
if ( strcmp(argv[4], "none") != 0 ) {
   rc = Coords_writeToFile(&coords, argv[4]) ;
   fprintf(msgFile, "\n return value %d from Coords_writeToFile(%p,%s)",
           rc, &coords, argv[4]) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
