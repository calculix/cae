/*  draw.c  */

#include "../Tree.h"

#define MYDEBUG 1

/*--------------------------------------------------------------------*/
static void
findLocalCoords ( int n, double x[], double xloc[], double rscale,
   double radius[], double xmin, double xmax ) ;
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to write an EPS file with a picture of a tree.
              each node can have its own radius and label

   filename  -- name of the file to be written
   xDV       -- x coordinates
   yDV       -- y coordinates
   rscale    -- scaling factor for radius of nodes
   radiusDV  -- radius of nodes, if NULL then radius = 1
   labelflag -- flag to specify whether labels are to be drawn
           1     -- draw labels
       otherwise -- do not draw labels
   fontscale -- scaling factor for font
   labelsIV  -- IV object that contains the labels of the nodes.
       if NULL then the node ids are used
   bbox[] -- bounding box for figure
      bbox[0] -- x_min
      bbox[1] -- y_min
      bbox[2] -- x_max
      bbox[3] -- y_max
   frame[] -- frame to hold tree
      frame[0] -- x_min
      frame[1] -- y_min
      frame[2] -- x_max
      frame[3] -- y_max
   bounds[] -- bounds for local coordinates
      if bounds is NULL then
         the tree fills the frame. note, this is a nonlinear process
         when the nodes have non-constant radii, and may not converge
         when the maximum radius is large when compared to the frame.
         if the process does not converge, a message is printed and
         the program exits.
      else
         bounds[0] -- xi_min
         bounds[1] -- eta_min
         bounds[2] -- xi_max
         bounds[3] -- eta_max
      endif

   recommendations, 
      bbox[] = { 0, 0, 500, 200 } for tall skinny trees
               { 0, 0, 500, 500 } for wide trees
      frame[0] = bbox[0] + 10
      frame[1] = bbox[1] + 10
      frame[2] = bbox[2] - 10
      frame[3] = bbox[3] - 10

   return value
      1 -- normal return
     -1 -- tree is NULL
     -2 -- filename is NULL
     -3 -- xDV is NULL
     -4 -- yDV is NULL
     -5 -- rscale is negative
     -6 -- fontscale is negative
     -7 -- bbox is NULL
     -8 -- frame is NULL

   created -- 99jan07, cca
   ------------------------------------------------------------
*/
int
Tree_drawToEPS (
   Tree     *tree,
   char     *filename,
   DV       *xDV,
   DV       *yDV,
   double   rscale,
   DV       *radiusDV,
   int      labelflag,
   double   fontscale,
   IV       *labelsIV,
   double   bbox[],
   double   frame[],
   double   bounds[]
) {
double   etamax, etamin, ximax, ximin, xmax, xmin, xrmax, xrmin,
         xscale, ymax, ymin, yrmax, yrmin, yscale ;
double   *radius, *x, *xloc, *y, *yloc ;
FILE     *fp ;
int      count, J, K, n ;
int      *fch, *par, *sib ;
/*
   ---------------
   check the input
   ---------------
*/
if ( tree == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n tree is NULL\n") ;
   return(-1) ;
}
if ( filename == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n filename is NULL\n") ;
   return(-2) ;
}
if ( xDV == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n xDV is NULL\n") ;
   return(-3) ;
}
if ( yDV == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n yDV is NULL\n") ;
   return(-4) ;
}
if ( rscale < 0.0 ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n rscale is negative\n") ;
   return(-5) ;
}
if ( fontscale < 0.0 ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n fontscale is negative\n") ;
   return(-6) ;
}
if ( bbox == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n bbox is NULL\n") ;
   return(-7) ;
}
if ( frame == NULL ) {
   fprintf(stderr, "\n error in Tree_drawToEPS()"
           "\n frame is NULL\n") ;
   return(-8) ;
}
n   = tree->n ;
par = tree->par ;
fch = tree->fch ;
sib = tree->sib ;
x   = DV_entries(xDV) ;
y   = DV_entries(yDV) ;
if ( radiusDV != NULL ) {
   radius = DV_entries(radiusDV) ;
} else {
   radius = NULL ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n\n x") ;
DVfprintf(stdout, n, x) ;
fprintf(stdout, "\n\n y") ;
DVfprintf(stdout, n, y) ;
if ( radius != NULL ) {
   fprintf(stdout, "\n\n radius") ;
   DVfprintf(stdout, n, radius) ;
}
#endif
xloc = DVinit(n, 0.0) ;
yloc = DVinit(n, 0.0) ;
if ( bounds != NULL ) {
/*
   ------------------------------------------
   get the local coordinates w.r.t the bounds
   ------------------------------------------
*/
   double   etamax, etamin, ximax, ximin, xmax, xmin, xoff, xscale,
            ymax, ymin, yoff, yscale ;
   xmin   = frame[0]  ; xmax   = frame[2]  ;
   ximin  = bounds[0] ; ximax  = bounds[2] ;
   xoff   = (xmin*ximax - xmax*ximin)/(ximax - ximin) ;
   xscale = (xmax - xmin)/(ximax - ximin) ;
   for ( J = 0 ; J < n ; J++ ) {
      xloc[J] = xoff + xscale*x[J] ;
   }
   ymin   = frame[1]  ; ymax   = frame[3]  ;
   etamin = bounds[1] ; etamax = bounds[3] ;
   yoff   = (ymin*etamax - ymax*etamin)/(etamax - etamin) ;
   yscale = (ymax - ymin)/(etamax - etamin) ;
   for ( J = 0 ; J < n ; J++ ) {
      yloc[J] = yoff + yscale*y[J] ;
   }
} else {
/*
   -----------------------------------------
   scale x[] and y[] to fit within the frame
   -----------------------------------------
*/
   xmin = frame[0] ;
   ymin = frame[1] ;
   xmax = frame[2] ;
   ymax = frame[3] ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n xmin = %.3g, xmax = %.3g", xmin, xmax) ;
#endif
   findLocalCoords(n, x, xloc, rscale, radius, xmin, xmax) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n\n ymin = %.3g, ymax = %.3g", ymin, ymax) ;
#endif
   findLocalCoords(n, y, yloc, rscale, radius, ymin, ymax) ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n xloc") ;
   DVfprintf(stdout, n, xloc) ;
#endif
#if MYDEBUG > 0
   fprintf(stdout, "\n\n yloc") ;
   DVfprintf(stdout, n, yloc) ;
#endif
/*
   -------------
   open the file
   -------------
*/
if ( (fp = fopen(filename, "w")) == NULL ) {
   fprintf(stderr, "\n unable to open file %s", filename) ;
   exit(-1) ;
}
/*
   ----------------------------
   print the header information
   ----------------------------
*/
fprintf(fp, 
        "%%!PS-Adobe-2.0 EPSF-1.2"
        "\n%%%%BoundingBox: %.3g %.3g %.3g %.3g",
        bbox[0], bbox[1], bbox[2], bbox[3]) ;
fprintf(fp, 
        "\n /CSH {"
        "\n %%"
        "\n %% center show a string"
        "\n %%"
        "\n %%  stack"
        "\n %%     string str"
        "\n %%"
        "\n dup stringwidth pop 2 div neg 0 rmoveto"
        "\n show"
        "\n } def") ;
fprintf(fp, 
        "\n /ML {"
        "\n %%"
        "\n %% moveto lineto"
        "\n %%"
        "\n %%  stack"
        "\n %%     x0 y0 x1 y1"
        "\n %%"
        "\n moveto lineto"
        "\n } def") ;
fprintf(fp, 
        "\n /FC {"
        "\n %%"
        "\n %% draw filled circle"
        "\n %%"
        "\n %%  stack"
        "\n %%     x y r"
        "\n %%"
        "\n newpath 2 index 1 index add 2 index moveto 0 360 arc fill"
        "\n } def") ;
fprintf(fp, 
        "\n /OC {"
        "\n %%"
        "\n %% draw open circle"
        "\n %%"
        "\n %%  stack"
        "\n %%     x y r"
        "\n %%"
        "\n newpath 2 index 1 index add 2 index moveto 0 360 arc stroke"
        "\n } def") ;
fprintf(fp, "\n /rscale    %.3f def", rscale) ;
fprintf(fp, "\n /fontscale %.3f def", fontscale) ;
/*
   --------------
   draw the edges
   --------------
*/
fprintf(fp, "\n %.3g %.3g %.3g %.3g rectclip",
        frame[0], frame[1], frame[2] - frame[0], frame[3] - frame[1]) ;
par = tree->par ;
count = 0 ;
for ( J = 0 ; J < n ; J++ ) {
   if ( (K = par[J]) != -1 ) {
      if ( count == 0 ) {
         fprintf(fp, "\n newpath") ;
      }
      fprintf(fp, "\n   %.3g %.3g %.3g %.3g ML",
              xloc[J], yloc[J], xloc[K], yloc[K]) ;
      count++ ;
      if ( count == 100 ) {
         fprintf(fp, "\n stroke") ;
         count = 0 ;
      }
   }
}
if ( count > 0 ) {
   fprintf(fp, "\n stroke") ;
}
/*
   -------------------------
   draw the nodes and labels
   -------------------------
*/
fprintf(fp, "\n\n gsave") ;
if ( labelflag == 1 ) {
   fprintf(fp, 
           "\n  /Helvetica-Bold findfont fontscale scalefont setfont") ;
}
if ( radius == NULL ) {
   for ( J = 0 ; J < n ; J++ ) {
      fprintf(fp, "\n    1.0 setgray") ;
      fprintf(fp, " %.3g %.3g %.3g FC", 
              xloc[J], yloc[J], rscale) ;
      fprintf(fp, "\n    0.0 setgray") ;
      fprintf(fp, " %.3g %.3g %.3g OC", 
              xloc[J], yloc[J], rscale) ;
      if ( labelflag == 1 ) {
         fprintf(fp, "\n   %.3g %.3g moveto ", 
                 xloc[J], yloc[J] - 0.5*rscale) ;
         if ( labelsIV != NULL ) {
            fprintf(fp, " (%d) CSH", IV_entry(labelsIV, J)) ;
         } else {
            fprintf(fp, " (%d) CSH", J) ;
         }
      }
   }
} else {
   for ( J = 0 ; J < n ; J++ ) {
      fprintf(fp, "\n    1.0 setgray") ;
      fprintf(fp, " %.3g %.3g %.3g FC", 
              xloc[J], yloc[J], rscale*radius[J]) ;
      fprintf(fp, "\n    0.0 setgray") ;
      fprintf(fp, " %.3g %.3g %.3g OC", 
              xloc[J], yloc[J], rscale*radius[J]) ;
      if ( labelflag == 1 ) {
         fprintf(fp, "\n   %.3g %.3g %.3g sub moveto ", 
                 xloc[J], yloc[J], 0.25*fontscale) ;
         if ( labelsIV != NULL ) {
            fprintf(fp, " (%d) CSH", IV_entry(labelsIV, J)) ;
         } else {
            fprintf(fp, " (%d) CSH", J) ;
         }
      }
   }
}
fprintf(fp, "\n\n grestore") ;
fprintf(fp, "\n %.3g %.3g %.3g %.3g rectstroke",
        frame[0], frame[1], frame[2] - frame[0], frame[3] - frame[1]) ;
fprintf(fp, "\n\n showpage") ;

return(1) ; }

/*--------------------------------------------------------------------*/
static void
findLocalCoords (
   int      n,
   double   x[],
   double   xloc[],
   double   rscale,
   double   radius[],
   double   xmin,
   double   xmax
) {
int      J, Jmax, Jmin ;
double   b1, b2, locwidth, width, ximax, ximin, 
         xlocmax, xlocmin, xoff, xscale ;

width = xmax - xmin ;
ximin = DVmin(n, x, &J) ;
ximax = DVmax(n, x, &J) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n ximin %f, ximax %f", ximin, ximax) ;
#endif
if ( ximax == ximin ) {
   DVfill(n, xloc, 0.5*(xmax+xmin)) ;
   return ;
}
xscale = width/(ximax - ximin) ;
#if MYDEBUG > 0
   fprintf(stdout, "\n initial xscale %f", xscale) ;
#endif
for ( J = 0 ; J < n ; J++ ) {
   xloc[J] = xmin + xscale*(x[J] - ximin) ;
}
while ( 1 ) {
   if ( radius == NULL ) {
      xlocmin = xloc[0] - rscale ;
      xlocmax = xloc[0] + rscale ;
      Jmin = Jmax = 0 ;
      for ( J = 0 ; J < n ; J++ ) {
         if ( xlocmin > xloc[J] - rscale ) {
              xlocmin = xloc[J] - rscale ;
              Jmin = J ;
         }
         if ( xlocmax < xloc[J] + rscale ) {
              xlocmax = xloc[J] + rscale ;
              Jmax = J ;
         }
      }
   } else {
      xlocmin = xloc[0] - rscale*radius[0] ;
      xlocmax = xloc[0] + rscale*radius[0] ;
      Jmin = Jmax = 0 ;
      for ( J = 0 ; J < n ; J++ ) {
         if ( xlocmin > xloc[J] - rscale*radius[J] ) {
              xlocmin = xloc[J] - rscale*radius[J] ;
              Jmin = J ;
         }
         if ( xlocmax < xloc[J] + rscale*radius[J] ) {
              xlocmax = xloc[J] + rscale*radius[J] ;
              Jmax = J ;
         }
      }
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n Jmin = %d, Jmax = %d", Jmin, Jmax) ;
   fprintf(stdout, "\n xlocmin %f, xlocmax %f", xlocmin, xlocmax) ;
#endif
   if ( Jmin == Jmax ) {
      DVfill(n, xloc, (xmin + xmax)/2) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n leaving") ;
#endif
      break ;
   } else {
      locwidth = xlocmax - xlocmin ;
#if MYDEBUG > 0
      fprintf(stdout, "\n width %f, locwidth %f", width, locwidth) ;
#endif
      if ( locwidth > 1.01*width || locwidth < 0.99*width ) {
         if ( radius == NULL ) {
            b1 = xmin + rscale ;
            b2 = xmax - rscale ;
         } else {
            b1 = xmin + rscale*radius[Jmin] ;
            b2 = xmax - rscale*radius[Jmax] ;
         }
         if ( b1 > b2 ) {
            fprintf(stderr, "\n\n error in Tree_drawEPS()"
                    "\n nonlinear process is unable to converge"
                    "\n reduce radius scaling factor\n") ;
            exit(-1) ;
         }
#if MYDEBUG > 0
         fprintf(stdout, "\n 1. x[%d] = %f, x[%d] = %f",
                 Jmin, x[Jmin], Jmax, x[Jmax]) ;
         fprintf(stdout, "\n 1. b1 = %f, b2 = %f", b1, b2) ;
#endif
         xscale = (b2 - b1)/(x[Jmax] - x[Jmin]) ;
         xoff   = -(b2*x[Jmin] - b1*x[Jmax])/(x[Jmax] - x[Jmin]) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n 1. xscale = %f, xoff = %f", xscale, xoff) ;
#endif
         for ( J = 0 ; J < n ; J++ ) {
            xloc[J] = xoff + xscale*x[J] ;
         }
      } else {
         break ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
