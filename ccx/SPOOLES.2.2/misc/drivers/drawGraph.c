/*  drawGraph.c  */

#include "../misc.h"
#include "../../Graph.h"
#include "../../Coords.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------
   draw a graph.

   (1) read a Graph object
   (2) read a Coords object
   (3) read an IV object that contains a tag vector.
       if (v,w) is an edge in the graph
       and tags[v] == tags[w] then
          draw edge (v,w)
   (4) bbox[4] is the bounding box for the plot.
       bbox = { xsw, ysw, xne, yne }
       try bbox = { 0, 0, 500, 500 }
       because coordinates are measured in points,
       72 points per inch.
   (5) rect[4] contains the frame for the plot.
       to put a 20 point margin around the plot,
       rect[0] = bbox[0] + 20
       rect[1] = bbox[1] + 20
       rect[2] = bbox[2] - 20
       rect[3] = bbox[3] - 20

   created -- 96jun06, cca
   -------------------------------------------------
*/
{
char     *fnCoords, *fnEPS, *fnGraph, *fnTagsIV ;
Coords   *coords ;
double   linewidth1, linewidth2, radius, t1, t2 ;
double   bbox[4], rect[4] ;
int      msglvl, nvtx, rc ;
IV       *tagsIV ;
Graph    *graph ;
FILE     *msgFile ;

if ( argc != 18 ) {
   fprintf(stdout, 
"\n\n usage : drawGraph msglvl msgFile GraphFile CoordsFile"
"\n         tagsIVfile epsFile bbox[4] rect[4] radius"
"\n    msglvl     -- message level"
"\n    msgFile    -- message file"
"\n    GraphFile  -- input graph file, must be *.graphf or *.graphb"
"\n    CoordsFile -- input Coords file, must be *.coordsf or *.coordsb"
"\n    tagsIVfile -- input IV file, must be *.ivf or *.ivb"
"\n                  contains the component ids"
"\n    epsFile    -- output postscript file, must be *.eps"
"\n    linewidth1 -- line width for edges connecting vertices"
"\n                  in the same component"
"\n    linewidth2 -- line width for edges connecting vertices"
"\n                  in the same component"
"\n    bbox[4]    -- bounding box for drawing"
"\n    rect[4]    -- rectangle to contain graph"
"\n    radius     -- radius for vertex circle"
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
fnGraph    = argv[3] ;
fnCoords   = argv[4] ;
fnTagsIV   = argv[5] ;
fnEPS      = argv[6] ;
linewidth1 = atof(argv[7]) ;
linewidth2 = atof(argv[8]) ;
bbox[0]    = atof(argv[9]) ;
bbox[1]    = atof(argv[10]) ;
bbox[2]    = atof(argv[11]) ;
bbox[3]    = atof(argv[12]) ;
rect[0]    = atof(argv[13]) ;
rect[1]    = atof(argv[14]) ;
rect[2]    = atof(argv[15]) ;
rect[3]    = atof(argv[16]) ;
radius     = atof(argv[17]) ;
fprintf(msgFile, 
        "\n drawGraph "
        "\n msglvl     -- %d" 
        "\n msgFile    -- %s" 
        "\n GraphFile  -- %s" 
        "\n CoordsFile -- %s" 
        "\n tagsIVfile -- %s" 
        "\n epsFile    -- %s" 
        "\n linewidth1 -- %f"
        "\n linewidth2 -- %f"
        "\n bbox[4] = { %f, %f, %f, %f }"
        "\n rect[4] = { %f, %f, %f, %f }"
        "\n radius     -- %f" 
        "\n",
        msglvl, argv[2], fnGraph, fnCoords, fnTagsIV, fnEPS,
        linewidth1, linewidth2, bbox[0], bbox[1], bbox[2], bbox[3],
        rect[0], rect[1], rect[2], rect[3], radius) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
graph = Graph_new() ;
if ( strcmp(fnGraph, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, fnGraph) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, fnGraph) ;
nvtx = graph->nvtx ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, fnGraph) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Graph object from file %s",
           argv[3]) ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the Coords object
   -------------------------
*/
if ( strcmp(fnCoords, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
coords = Coords_new() ;
rc = Coords_readFromFile(coords, fnCoords) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in coords from file %s",
        t2 - t1, fnCoords) ;
if ( rc != 1 ) {
  fprintf(msgFile, "\n return value %d from Coords_readFromFile(%p,%s)",
           rc, coords, fnCoords) ;
   exit(-1) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n after reading Coords object from file %s",
           argv[3]) ;
   Coords_writeForHumanEye(coords, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   read in the IV object
   ---------------------
*/
if ( strcmp(fnTagsIV, "none") == 0 ) {
   tagsIV = NULL ;
} else {
   MARKTIME(t1) ;
   tagsIV = IV_new() ;
   rc = IV_readFromFile(tagsIV, fnTagsIV) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : read in tagsIV from file %s",
           t2 - t1, fnTagsIV) ;
   if ( rc != 1 ) {
     fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
              rc, tagsIV, fnTagsIV) ;
      exit(-1) ;
   }
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n\n after reading IV object from file %s",
              argv[3]) ;
      IV_writeForHumanEye(tagsIV, msgFile) ;
      fflush(msgFile) ;
   }
}
fprintf(msgFile, "\n getting ready to call drawGraphEPS") ;
fflush(msgFile) ;
/*
   --------------
   draw the graph
   --------------
*/
drawGraphEPS(graph, coords, tagsIV, bbox, rect, linewidth1, linewidth2, 
             radius, fnEPS, msglvl, msgFile) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;
Coords_free(coords) ;
if ( tagsIV != NULL ) {
   IV_free(tagsIV) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
