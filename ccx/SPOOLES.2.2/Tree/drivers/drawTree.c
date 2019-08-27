/*  drawTree.c  */

#include "../Tree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------
   draw the tree

   created -- 99jan23, cca
   ----------------------------------------
*/
{
char     coordflag, heightflag ;
char     *inTagsFileName, *inTreeFileName, *outEPSfileName ;
double   fontsize, radius, t1, t2 ;
double   bbox[4], frame[4] ;
DV       *xDV, *yDV ;
int      ierr, msglvl, rc, tagsflag ;
IV       *tagsIV ;
Tree     *tree ;
FILE     *msgFile ;

if ( argc != 19 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inTreeFile inTagsFile outEPSfile "
"\n       heightflag coordflag radius bbox[4] frame[4] tagflag fontsize"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inTreeFile -- input file, must be *.treef or *.treeb"
      "\n    inTagsFile -- input file, must be *.ivf or *.ivb or none"
      "\n    outEPSfile -- output file"
      "\n    heightflag -- height flag"
      "\n       'D' -- use depth metric"
      "\n       'H' -- use height metric"
      "\n    coordflag -- coordinate flag"
      "\n       'C' -- use (x,y) Cartesian coordinates"
      "\n       'P' -- use (r,theta) polar coordinates"
      "\n    radius   -- radius of node"
      "\n    bbox[4]  -- bounding box"
      "\n    frame[4] -- frame for plot"
      "\n    fontsize -- size of fonts (in points)"
      "\n    tagflag  -- if 1, draw labels, otherwise, do not"
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
inTreeFileName = argv[3] ;
inTagsFileName = argv[4] ;
outEPSfileName = argv[5] ;
heightflag     = argv[6][0] ;
coordflag      = argv[7][0] ;
radius         = atof(argv[8]) ;
bbox[0]        = atof(argv[9]) ;
bbox[1]        = atof(argv[10]) ;
bbox[2]        = atof(argv[11]) ;
bbox[3]        = atof(argv[12]) ;
frame[0]       = atof(argv[13]) ;
frame[1]       = atof(argv[14]) ;
frame[2]       = atof(argv[15]) ;
frame[3]       = atof(argv[16]) ;
fontsize       = atof(argv[17]) ;
tagsflag       = atoi(argv[18]) ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl     -- %d" 
        "\n msgFile    -- %s" 
        "\n inTreeFile -- %s" 
        "\n inTagsFile -- %s" 
        "\n outEPSfile -- %s" 
        "\n heightflag -- %c" 
        "\n coordflag  -- %d" 
        "\n radius     -- %.3g" 
        "\n bbox       -- %.3g %.3g %.3g %.3g" 
        "\n frame      -- %.3g %.3g %.3g %.3g" 
        "\n fontsize   -- %.3g"
        "\n",
        argv[0], msglvl, argv[2], inTreeFileName, inTagsFileName,
        outEPSfileName, heightflag, coordflag, radius, 
        bbox[0], bbox[1], bbox[2], bbox[3],
        frame[0], frame[1], frame[2], frame[3], fontsize, tagsflag) ;
fflush(msgFile) ;
/*
   ------------------------
   read in the Tree object
   ------------------------
*/
if ( strcmp(inTreeFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
tree = Tree_new() ;
MARKTIME(t1) ;
rc = Tree_readFromFile(tree, inTreeFileName) ;
/*
Tree_setFchSibRoot(tree) ;
*/
Tree_leftJustify(tree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in tree from file %s",
        t2 - t1, inTreeFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from Tree_readFromFile(%p,%s)",
           rc, tree, inTreeFileName) ;
   exit(-1) ;
}
fprintf(msgFile, "\n\n after reading Tree object from file %s",
        inTreeFileName) ;
if ( msglvl > 2 ) {
   Tree_writeForHumanEye(tree, msgFile) ;
} else {
   Tree_writeStats(tree, msgFile) ;
}
fflush(msgFile) ;
if ( Tree_maxNchild(tree) > 2 ) {
   fprintf(msgFile, "\n\n maximum number of children = %d",
           Tree_maxNchild(tree)) ;
}
if ( strcmp(inTagsFileName, "none") != 0 ) {
/*
   --------------------------
   read in the tags IV object
   --------------------------
*/
   tagsIV = IV_new() ;
   MARKTIME(t1) ;
   rc = IV_readFromFile(tagsIV, inTagsFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : read in tagsIV from file %s",
           t2 - t1, inTagsFileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
              rc, tagsIV, inTagsFileName) ;
      exit(-1) ;
   }
   fprintf(msgFile, "\n\n after reading IV object from file %s",
           inTagsFileName) ;
   if ( msglvl > 2 ) {
      IV_writeForHumanEye(tagsIV, msgFile) ;
   } else {
      IV_writeStats(tagsIV, msgFile) ;
   }
   fflush(msgFile) ;
   if ( IV_size(tagsIV) != tree->n ) {
      fprintf(stderr, 
              "\n fatal error, IV_size(tagsIV) = %d, tree->n = %d",
              IV_size(tagsIV), tree->n) ;
      exit(-1) ;
   }
} else {
   tagsIV = NULL ;
}
/*
   -------------------------------
   get the coordinates of the tree
   -------------------------------
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
   fflush(msgFile) ;
}
/*
   -------------
   draw the Tree
   -------------
*/
rc = Tree_drawToEPS(tree, outEPSfileName, xDV, yDV, radius, NULL,
                    tagsflag, fontsize, tagsIV, bbox, frame, NULL) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from Tree_drawToEPSfile()", rc) ;
   exit(-1) ;
}
/*
   ---------------------
   free the Tree object
   ---------------------
*/
Tree_free(tree) ;
if ( tagsIV != NULL ) {
   IV_free(tagsIV) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
