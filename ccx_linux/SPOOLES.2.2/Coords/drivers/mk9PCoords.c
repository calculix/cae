/*  mk9PCoords.c  */

#include "../Coords.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/

int
main ( int argc, char *argv[] )
/*
   ------------------------------------
   create a Coords object for a 9P grid
   
   created -- 96feb02, cca
   ------------------------------------
*/
{
Coords   *coords ;
double   t1, t2 ;
FILE     *msgFile ;
float    bbox[4] = { 0.0, 0.0, 1.0, 1.0 } ;
int      msglvl, n1, n2 ;

if ( argc != 6 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile n1 n2 outCoordsFile \n"
"\n msglvl        -- message level"
"\n msgFile       -- message file"
"\n n1            -- # of grid points in first direction"
"\n n2            -- # of grid points in second direction"
"\n outCoordsFile -- file to write out Coords object"
"\n", argv[0]) ;
   return(1) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s", argv[0], argv[2]) ;
   exit(-1) ; 
}
n1 = atoi(argv[3]) ;
n2 = atoi(argv[4]) ;
fprintf(msgFile, "\n\n %s %d %s %d %d %s",
        argv[0], msglvl, argv[2], n1, n2, argv[5]) ;
fflush(msgFile) ;
/*
   ------------------------
   create the Coords object
   ------------------------
*/
MARKTIME(t1) ;
coords = Coords_new() ;
Coords_init9P(coords, bbox, COORDS_BY_TUPLE, n1, n2, 1) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %.3f : create Coords object", (t2 - t1)) ;
if ( msglvl <= 2 ) {
   Coords_writeStats(coords, msgFile) ;
   fflush(msgFile) ;
} else {
   Coords_writeForHumanEye(coords, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------------
   write out the Coords object
   ---------------------------
*/
if ( strcmp("none", argv[5]) != 0 ) {
   MARKTIME(t1) ;
   Coords_writeToFile(coords, argv[5]) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %.3f : write out Coords to file %s",
           (t2 - t1), argv[5]) ;
}

fprintf(msgFile, "\n") ;
fflush(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
