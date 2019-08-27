/*  drawGraphEPS.c  */

#include "../../Graph.h"
#include "../../Coords.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   draw a graph to an EPS file

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

   created -- 98apr11, cca
   -------------------------------------------------
*/
void
drawGraphEPS (
   Graph    *graph,
   Coords   *coords,
   IV       *tagsIV,
   double   bbox[],
   double   rect[],
   double   linewidth1,
   double   linewidth2,
   double   radius,
   char     *epsFileName,
   int      msglvl,
   FILE     *msgFile
) {
double   a, b, d, height, offset, width, 
         xmax, xmin, xsize, xv, xw, x0, x1, 
         ymax, ymin, ysize, yv, yw, y0, y1 ;
FILE     *epsFile ;
int      ii, nedge, nvtx, v, vsize, w ;
int      *tags, *vadj ;

nvtx = graph->nvtx ;
if ( tagsIV == NULL ) {
   tags = NULL ;
} else {
   tags = IV_entries(tagsIV) ;
}
/*
   -----------------
   open the EPS file
   -----------------
*/
if ( strcmp(epsFileName, "stdout") == 0 ) {
   epsFile = stdout ;
} else if ( (epsFile = fopen(epsFileName, "w")) == NULL ) {
   fprintf(stderr, "\n fatal error in drawGraphEPS"
           "\n unable to open file %s\n", epsFileName) ;
   return ;
}
/*
   -----------------------------------
   write the preamble for the EPS file
   -----------------------------------
*/
fprintf(epsFile,
   "%%!PS-Adobe-2.0 EPSF-1.2"
   "\n%%%%BoundingBox: %.1f %.1f %.1f %.1f",
   bbox[0], bbox[1], bbox[2], bbox[3]) ;

fprintf(epsFile,
   "\n /radius %.3f def"
   "\n /Helvetica findfont %.3f scalefont setfont"
   "\n /M {moveto} def"
   "\n /L {lineto} def"
   "\n /ACF { %% stack : x y radius"
   "\n    newpath 0 360 arc closepath fill "
   "\n } def"
   "\n /str 6 string def"
   "\n /drawLabel { %% x y label radius"
   "\n    /radius exch def"
   "\n    /label  exch def"
   "\n    /y      exch def"
   "\n    /x      exch def"
   "\n    gsave"
   "\n       1.0 setgray"
   "\n       x radius add y moveto"
   "\n       x y radius 0 360 arc"
   "\n       fill"
   "\n       0.0 setgray"
   "\n       x radius add y moveto"
   "\n       x y radius 0 360 arc"
   "\n       stroke"
   "\n       x y moveto"
   "\n       label stringwidth pop 2 div neg radius 2 div neg rmoveto"
   "\n       label show"
   "\n    grestore"
   "\n } def ", radius, 1.25*radius) ;
/*
   ---------------------------------------
   determine the transformation parameters
   ---------------------------------------
*/
xmin   = Coords_min(coords, 1) ;
xmax   = Coords_max(coords, 1) ;
ymin   = Coords_min(coords, 2) ;
ymax   = Coords_max(coords, 2) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
           "\n xmin = %.3g, xmax = %.3g, ymin = %.3g, ymax = %.3g", 
           xmin, xmax, ymin, ymax) ;
}
xsize  = xmax - xmin ;
ysize  = ymax - ymin ;
width  = rect[2] - rect[0] ;
height = rect[3] - rect[1] ;
if ( msglvl > 2 ) {
   fprintf(msgFile, 
        "\n xsize = %.3g, ysize = %.3g, width = %.3g, height = %.3g", 
        xsize, ysize, width, height) ;
}
if ( ysize * width <= xsize * height ) {
   a = width / xsize ;
   b = rect[0] ;
   offset = (rect[3] - rect[1] - a * ysize)/2 ;
   d = rect[1] + offset - a * ymin ;
} else {
   a = height / ysize ;
   d = rect[1] ;
   offset = (rect[2] - rect[0] - a * xsize)/2 ;
   b = rect[0] + offset - a * xmin ;
}
if ( ysize * width <= xsize * height ) {
   a = width / xsize ;
} else {
   a = height / ysize ;
}
b = 0.5*(rect[2] + rect[0] - a*(xmin + xmax)) ;
d = 0.5*(rect[3] + rect[1] - a*(ymin + ymax)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n width = %.3g, height = %.3g", width, height) ;
   fprintf(msgFile, "\n xsize = %.3g, ysize = %.3g", xsize, ysize) ;
   fprintf(msgFile, 
           "\n xmin = %.3g, xmax = %.3g, ymin = %.3g, ymax = %.3g",
           xmin, xmax, ymin, ymax) ;
   fprintf(msgFile, "\n a = %.3g, b = %.3g, d = %.3g", a, b, d) ;
}
if ( tags == NULL ) {
/*
   --------------------------------
   no component ids, draw the edges
   --------------------------------
*/
   fprintf(epsFile,
           "\n gsave"
           "\n   %.3f setlinewidth"
           "\n   0.0 setgray", linewidth1) ;
   nedge = 0 ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      Graph_adjAndSize(graph, v, &vsize, &vadj) ;
      xv = Coords_value(coords, 1, v) ;
      yv = Coords_value(coords, 2, v) ;
      x0 = a * xv + b ;
      y0 = a * yv + d ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         if ( w < v ) {
            xw = Coords_value(coords, 1, w) ;
            yw = Coords_value(coords, 2, w) ;
            x1 = a * xw + b ;
            y1 = a * yw + d ;
            if ( nedge % 100 == 0 ) {
               fprintf(epsFile, "\n    newpath") ;
            }
            fprintf(epsFile, "\n       %.3g %.3g M %.3g %.3g L", 
                    x0, y0, x1, y1) ;
            if ( ++nedge % 100 == 0 ) {
               fprintf(epsFile, "\n    stroke") ;
            }
         }
      }
   }
   if ( nedge % 100 != 0 ) {
      fprintf(epsFile, "\n    stroke") ;
   }
   fprintf(epsFile,
           "\n grestore") ;
   fprintf(epsFile,
           "\n gsave"
           "\n   0.1 setlinewidth"
           "\n   0.0 setgray") ;
   if ( radius > 0.0 ) {
/*
      -----------------
      draw the vertices
      -----------------
*/
      for ( v = 0 ; v < nvtx ; v++ ) {
         xv = Coords_value(coords, 1, v) ;
         yv = Coords_value(coords, 2, v) ;
         x0 = a * xv + b ;
         y0 = a * yv + d ;
         fprintf(epsFile, "\n %.3f %.3f () radius drawLabel", 
                 x0, y0) ;
      }
   }
   fprintf(epsFile, "\n grestore") ;
} else {
/*
   -----------------------------------------
   component ids are present, draw the edges 
   between vertices in the same component
   -----------------------------------------
*/
   fprintf(epsFile,
           "\n gsave"
           "\n   %.3f setlinewidth"
           "\n   0.0 setgray", linewidth1) ;
   nedge = 0 ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      if ( tags[v] >= 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         xv = Coords_value(coords, 1, v) ;
         yv = Coords_value(coords, 2, v) ;
         x0 = a * xv + b ;
         y0 = a * yv + d ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( w < v && tags[w] == tags[v] ) {
               xw = Coords_value(coords, 1, w) ;
               yw = Coords_value(coords, 2, w) ;
               x1 = a * xw + b ;
               y1 = a * yw + d ;
               if ( nedge % 100 == 0 ) {
                  fprintf(epsFile, "\n    newpath") ;
               }
               fprintf(epsFile, "\n       %.3g %.3g M %.3g %.3g L", 
                       x0, y0, x1, y1) ;
               if ( ++nedge % 100 == 0 ) {
                  fprintf(epsFile, "\n    stroke") ;
               }
            }
         }
      }
   }
   if ( nedge % 100 != 0 ) {
      fprintf(epsFile, "\n    stroke") ;
   }
   fprintf(epsFile,
        "\n grestore") ;
   fprintf(epsFile,
           "\n gsave"
           "\n   %.3f setlinewidth"
           "\n   0.0 setgray", linewidth2) ;
   nedge = 0 ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      if ( tags[v] >= 0 ) {
         Graph_adjAndSize(graph, v, &vsize, &vadj) ;
         xv = Coords_value(coords, 1, v) ;
         yv = Coords_value(coords, 2, v) ;
         x0 = a * xv + b ;
         y0 = a * yv + d ;
         for ( ii = 0 ; ii < vsize ; ii++ ) {
            w = vadj[ii] ;
            if ( w < v && tags[w] != tags[v] && tags[w] >= 0 ) {
               xw = Coords_value(coords, 1, w) ;
               yw = Coords_value(coords, 2, w) ;
               x1 = a * xw + b ;
               y1 = a * yw + d ;
               if ( nedge % 100 == 0 ) {
                  fprintf(epsFile, "\n    newpath") ;
               }
               fprintf(epsFile, "\n       %.3g %.3g M %.3g %.3g L", 
                       x0, y0, x1, y1) ;
               if ( ++nedge % 100 == 0 ) {
                  fprintf(epsFile, "\n    stroke") ;
               }
            }
         }
      }
   }
   if ( nedge % 100 != 0 ) {
      fprintf(epsFile, "\n    stroke") ;
   }
   fprintf(epsFile,
        "\n grestore") ;
   fprintf(epsFile,
           "\n gsave"
           "\n   0.1 setlinewidth"
           "\n   0.0 setgray") ;
   if ( radius > 0.0 ) {
/*
      -----------------
      draw the vertices
      -----------------
*/
      for ( v = 0 ; v < nvtx ; v++ ) {
         if ( tags[v] >= 0 ) {
            xv = Coords_value(coords, 1, v) ;
            yv = Coords_value(coords, 2, v) ;
            x0 = a * xv + b ;
            y0 = a * yv + d ;
            fprintf(epsFile, "\n %.3f %.3f (%d) radius drawLabel", 
                    x0, y0, tags[v]) ;
         }
      }
   }
   fprintf(epsFile, "\n grestore") ;
}
fprintf(epsFile, "\n showpage") ;
/*
   ----------------------------
   close the file if not stdout
   ----------------------------
*/
if ( strcmp(epsFileName, "stdout") != 0 ) {
   fclose(epsFile) ;
}

return ; }

/*--------------------------------------------------------------------*/
