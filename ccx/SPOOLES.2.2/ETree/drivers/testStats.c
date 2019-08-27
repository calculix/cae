/*  testStats.c  */

#include "../../SymbFac.h"
#include "../../ETree.h"
#include "../../InpMtx.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------
   (1) read in an ETree object.
   (2) read in a Graph object.
   (3) draw a tree picture of the distribution
       of some metric across the fronts
   
   created -- 99jan22, cca
   -------------------------------------------
*/
{
char     coordflag, heightflag ;
char     *inETreeFileName, *inGraphFileName, *outEPSfileName ;
double   fontscale, rmax, rscale, t1, t2 ;
double   bbox[4], frame[4] ;
double   *radius ;
DV       *radiusDV, *xDV, *yDV ;
Graph    *graph ;
int      J, labelflag, metricType, msglvl, nfront, nvtx, rc ;
int      *metric, *vwghts ;
IV       *metricIV ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 12 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inGraphFile outEPSfile "
"\n         metricType heightflag coordflag rscale labelflag fontscale"
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
"\n    outEPSfile   -- output file for tree picture, must be *.eps"
"\n    metricType   -- metric type"
"\n       0 -- no metric"
"\n       1 -- # of nodes in front"
"\n       2 -- # of original matrix entries in front"
"\n       3 -- # of factor matrix entries in front"
"\n       4 -- # of forward factor ops in front"
"\n       5 -- # of backward factor ops in front"
"\n    heightflag   -- flag for type of height"
"\n      'D' -- nodes placed by depth in tree"
"\n      'H' -- nodes placed by height in tree"
"\n    coordflag    -- flag for type of coordinate"
"\n      'C' -- cartesian (x,y) placement"
"\n      'P' -- polar (r,theta) placement"
"\n    rmax      -- maximum radius in pts"
"\n    labelflag -- flag for type of labels"
"\n           1     -- draw labels"
"\n       otherwise -- do not draw labels"
"\n    fontscale -- scaling parameter for fonts when labels drawed"
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
inETreeFileName  = argv[3] ;
inGraphFileName  = argv[4] ;
outEPSfileName   = argv[5] ;
metricType       = atoi(argv[6]) ;
heightflag       = *argv[7] ;
coordflag        = *argv[8] ;
rmax             = atof(argv[9]) ;
labelflag        = atoi(argv[10]) ;
fontscale        = atof(argv[11]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inGraphFile   -- %s" 
        "\n outEPSfile    -- %s" 
        "\n metricType    -- %d" 
        "\n heightflag    -- %c" 
        "\n coordflag     -- %c" 
        "\n rmax          -- %f" 
        "\n labelflag     -- %d" 
        "\n fontscale     -- %f" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, inGraphFileName, 
        outEPSfileName, metricType, heightflag, coordflag, rmax, 
        labelflag, fontscale) ;
fflush(msgFile) ;
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
ETree_leftJustify(etree) ;
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
fflush(msgFile) ;
nfront = ETree_nfront(etree) ;
tree   = ETree_tree(etree) ;
/*
   ------------------------
   read in the Graph object
   ------------------------
*/
if ( strcmp(inGraphFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
graph = Graph_new() ;
MARKTIME(t1) ;
rc = Graph_readFromFile(graph, inGraphFileName) ;
nvtx = graph->nvtx ;
vwghts = graph->vwghts ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in graph from file %s",
        t2 - t1, inGraphFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Graph_readFromFile(%p,%s)",
           rc, graph, inGraphFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Graph object from file %s",
        inGraphFileName) ;
if ( msglvl > 2 ) {
   Graph_writeForHumanEye(graph, msgFile) ;
} else {
   Graph_writeStats(graph, msgFile) ;
}
fflush(msgFile) ;
/*
   ----------------------
   compute the statistics
   ----------------------
*/
metricIV = IV_new() ;
IV_init(metricIV, nfront, NULL) ;
metric = IV_entries(metricIV) ;
switch ( metricType ) {
case 0 : {
/*
   ---------
   no metric
   ---------
*/
   IVfill(nfront, metric, 1) ;
   } break ;
case 1 : {
/*
   -----------
   front sizes
   -----------
*/
   IVcopy(nfront, metric, ETree_nodwghts(etree)) ;
   } break ;
case 2 : {
/*
   ----------------------------------------------
   number of original matrix entries in the front
   ----------------------------------------------
*/
   int   I, ii, J, v, vsize, w ;
   int   *vadj, *vtxToFront ;

   vtxToFront = ETree_vtxToFront(etree) ;
   for ( v = 0 ; v < nvtx ; v++ ) {
      I = vtxToFront[v] ;
      Graph_adjAndSize(graph, v, &vsize, &vadj) ;
      for ( ii = 0 ; ii < vsize ; ii++ ) {
         w = vadj[ii] ;
         J = vtxToFront[w] ;
         if ( I <= J ) {
            if ( vwghts == NULL ) {
               metric[I]++ ;
            } else {
               metric[I] += vwghts[v]*vwghts[w] ;
            }
         } else {
            if ( vwghts == NULL ) {
               metric[J]++ ;
            } else {
               metric[J] += vwghts[v]*vwghts[w] ;
            }
         }
      }
   }
   } break ;
case 3 : {
/*
   --------------------------------------------
   number of factor matrix entries in the front
   --------------------------------------------
*/
   int   *nzf ;
   IV    *nzfIV ;

   nzfIV = ETree_factorEntriesIV(etree, SPOOLES_SYMMETRIC) ;
   nzf   = IV_entries(nzfIV) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      metric[J] = nzf[J] ;
   }
   IV_free(nzfIV) ;
   } break ;
case 4 : {
/*
   ------------------------------------------------
   number of forward factor operations in the front
   ------------------------------------------------
*/
   double   *ops ;
   DV       *opsDV ;

   opsDV = ETree_forwardOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
   ops   = DV_entries(opsDV) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      metric[J] = ops[J] ;
   }
   DV_free(opsDV) ;
   } break ;
case 5 : {
/*
   -------------------------------------------------
   number of backward factor operations in the front
   -------------------------------------------------
*/
   double   *ops ;
   DV       *opsDV ;
   IVL      *symbfacIVL ;

   symbfacIVL = SymbFac_initFromGraph(etree, graph) ;
   opsDV = ETree_backwardOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC,
                             graph->vwghts, symbfacIVL) ;
   ops   = DV_entries(opsDV) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      metric[J] = ops[J] ;
   }
   DV_free(opsDV) ;
   IVL_free(symbfacIVL) ;
   } break ;
default :
   fprintf(stderr, "\n error in testStats"
           "\n metricType %d not supported", metricType) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n metric vector") ;
   IV_writeForHumanEye(metricIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------------
   set the radius vector
   ---------------------
*/
radiusDV = DV_new() ;
DV_init(radiusDV, nfront, NULL) ;
radius = DV_entries(radiusDV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   radius[J] = sqrt(0.318309886*metric[J]) ;
}
rscale = rmax / DVmax(nfront, radius, &J) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n radius vector") ;
   DV_writeForHumanEye(radiusDV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ------------------------
   get the tree coordinates
   ------------------------
*/
xDV = DV_new() ;
yDV = DV_new() ;
rc = Tree_getSimpleCoords(tree, heightflag, coordflag, xDV, yDV) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from Tree_getSimpleCoords()",rc);
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n x-coordinates") ;
   DV_writeForHumanEye(xDV, msgFile) ;
   fprintf(msgFile, "\n\n y-coordinates") ;
   DV_writeForHumanEye(yDV, msgFile) ;
   fprintf(msgFile, "\n\n radius") ;
   DV_writeForHumanEye(radiusDV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------
   draw the tree
   -------------
*/
bbox[0]   =  50 ;
bbox[1]   =  50 ;
bbox[2]   = 550 ;
bbox[3]   = 550 ;
frame[0]  =  60 ;
frame[1]  =  60 ;
frame[2]  = 540 ;
frame[3]  = 540 ;
rc = Tree_drawToEPS(tree, outEPSfileName, xDV, yDV, rscale, radiusDV, 
                    labelflag, fontscale, NULL, bbox, frame, NULL) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from Tree_drawToEPSfile()", rc) ;
   exit(-1) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
Graph_free(graph) ;
IV_free(metricIV) ;
DV_free(radiusDV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
