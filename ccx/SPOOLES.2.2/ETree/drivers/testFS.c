/*  testFS.c  */

#include "../../ETree.h"
#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   read in an ETree object. compute the heights of the nodes 
   in the tree w.r.t. an out-of-core forward sparse factorization. 
   draw the tree.
   
   created -- 99jan07, cca
   -------------------------------------------------------
*/
{
char     *firstEPSfilename, *inETreeFileName, *secondEPSfilename ;
double   nfops1, radius, t1, t2 ;
double   *x, *y ;
double   bbox[4], bounds[4], frame[4] ;
DV       *xDV, *yDV ;
IV       *dmetricIV, *hmetricIV, *vmetricIV ;
int      bndJ, J, K, labelflag, maxdepth, maxnent, msglvl, 
         nfent1, nfind1, nfront, nleaves1, nnode1, rc ;
int      *bndwghts, *depths, *heights, *nzfs, *par ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 8 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile labelflag radius"
"\n         firstEPSfileName secondEPSfileName"
"\n    msglvl            -- message level"
"\n    msgFile           -- message file"
"\n    inETreeFile       -- input file, must be *.etreef or *.etreeb"
"\n    labelflag         -- flag to draw labels"
"\n    radius            -- radius of node"
"\n    firstEPSfilename  -- EPS file for subtree working storage"
"\n    secondEPSfilename -- EPS file for node working storage"
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
inETreeFileName   = argv[3] ;
labelflag         = atoi(argv[4]) ;
radius            = atof(argv[5]) ;
firstEPSfilename  = argv[6] ;
secondEPSfilename = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n firstEPSfilename    -- %s" 
        "\n labelflag     -- %d" 
        "\n radius        -- %f" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, firstEPSfilename,
        labelflag, radius) ;
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
fprintf(msgFile, "\n\n after reading ETree object from file %s",
        inETreeFileName) ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
ETree_leftJustify(etree) ;
fprintf(msgFile, "\n\n %d LU entries", ETree_nFactorEntries(etree, 2)) ;
fflush(msgFile) ;
/*
   ----------------------
   compute the statistics
   ----------------------
*/
tree   = etree->tree   ;
nfront = etree->nfront ;
nnode1 = etree->tree->n ;
nfind1 = ETree_nFactorIndices(etree) ;
nfent1 = ETree_nFactorEntries(etree, 1) ;
nfops1 = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
nleaves1 = Tree_nleaves(etree->tree) ;
fprintf(msgFile, "\n root front %d has %d vertices",
        etree->tree->root,
        etree->nodwghtsIV->vec[etree->tree->root]) ;
fprintf(msgFile, "\n %d fronts, %d indices, %d entries, %.0f ops",
        nfront, nfind1, nfent1, nfops1) ;
fprintf(msgFile, "\n max front size = %d",
        IV_max(ETree_nodwghtsIV(etree))) ;
fprintf(msgFile, "\n max boundary size = %d",
        IV_max(ETree_bndwghtsIV(etree))) ;
/*
   ------------------------------
   get the # of entries per front
   ------------------------------
*/
vmetricIV = ETree_factorEntriesIV(etree, SPOOLES_SYMMETRIC) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n entries per front") ;
   IV_writeForHumanEye(vmetricIV, msgFile) ;
   fflush(msgFile) ;
}
maxnent = IV_max(vmetricIV) ;
fprintf(msgFile, "\n\n max entries per front = %d", maxnent) ;
fflush(msgFile) ;
/*
   --------------------------------
   get the simple (x,y) coordinates
   --------------------------------
*/
xDV = DV_new() ;
yDV = DV_new() ;
rc = Tree_getSimpleCoords(tree, 'H', 'C', xDV, yDV) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in %s",
           "\n return value %d from Tree_getSimpleCoords()\n",
           argv[0], rc) ;
   exit(-1) ;
}
x = DV_entries(xDV) ;
y = DV_entries(yDV) ;
/*
   ----------------------
   get the height profile
   ----------------------
*/
hmetricIV = Tree_setHeightImetric(tree, vmetricIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n entries height per front") ;
   IV_writeForHumanEye(hmetricIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   -----------------------------------------------
   compute y[J] = heights[J] + |bndJ|*(|bndJ|+1)/2
   -----------------------------------------------
*/
heights = IV_entries(hmetricIV) ;
bndwghts = ETree_bndwghts(etree) ;
for ( J = 0 ; J < nfront ; J++ ) {
   bndJ = bndwghts[J] ;
   y[J] = heights[J] + (bndJ*(bndJ+1))/2 ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n  J    x    y") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      fprintf(msgFile, "\n %5d %12.3f %12.0f", J, x[J], y[J]) ;
   }
}
/*
   ------------------
   compute the bounds
   ------------------
*/
bounds[0] = 0.0 ;
bounds[1] = 0.0 ;
bounds[2] = DV_max(xDV) ;
bounds[3] = DV_max(yDV) ;
fprintf(stdout, "\n\n bounds = [ %.3g %.3g %.3g %.3g ] ",
        bounds[0], bounds[1], bounds[2], bounds[3]) ;
/*
   -------------
   draw the tree
   -------------
*/
bbox[0]  = 50.0 ; 
bbox[1]  = 50.0 ; 
bbox[2]  = 500.0 ; 
bbox[3]  = 500.0 ;
frame[0] = bbox[0] + 10 ;
frame[1] = bbox[1] + 10 ;
frame[2] = bbox[2] - 10 ;
frame[3] = bbox[3] - 10 ;
rc = Tree_drawToEPS(tree, firstEPSfilename, xDV, yDV, radius, NULL,
                    labelflag, radius, NULL, bbox, frame, bounds) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in %s"
           "\n return value %d from Tree_drawToEPS()\n",
           argv[0], rc) ;
   exit(-1) ;
}
/*
   ---------------------
   get the depth profile
   ---------------------
*/
dmetricIV = Tree_setDepthImetric(tree, vmetricIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n entries depth per front") ;
   IV_writeForHumanEye(dmetricIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   compute y[J] = depths[par(J)] + nfent[J]
   ----------------------------------------
*/
par    = Tree_par(tree) ;
depths = IV_entries(dmetricIV) ;
nzfs   = IV_entries(vmetricIV) ;
for ( J = 0 ; J < nfront ; J++ ) {
   y[J] = nzfs[J] ;
   if ( (K = par[J]) != -1 ) {
      y[J] += depths[K] ;
    }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n  J    x    y") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      fprintf(msgFile, "\n %5d %12.3f %12.0f", J, x[J], y[J]) ;
   }
}
/*
   -------------
   draw the tree
   -------------
*/
bbox[0]  = 50.0 ; 
bbox[1]  = 50.0 ; 
bbox[2]  = 500.0 ; 
bbox[3]  = 500.0 ;
frame[0] = bbox[0] + 10 ;
frame[1] = bbox[1] + 10 ;
frame[2] = bbox[2] - 10 ;
frame[3] = bbox[3] - 10 ;
rc = Tree_drawToEPS(tree, secondEPSfilename, xDV, yDV, radius, NULL,
                    labelflag, radius, NULL, bbox, frame, bounds) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error in %s"
           "\n return value %d from Tree_drawToEPSfile()\n",
           argv[0], rc) ;
   exit(-1) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
IV_free(vmetricIV) ;
IV_free(hmetricIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
