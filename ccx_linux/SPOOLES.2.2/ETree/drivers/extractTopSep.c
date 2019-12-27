/*  extractTopSep.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   read in a ETree object, create an IV object with the same size,
   mark the vertices in the top level separator(s), write the IV
   object to a file

   created -- 96may02, cca
   ---------------------------------------------------------------
*/
{
char     *inETreeFileName, *outIVfileName ;
double   t1, t2 ;
int      msglvl, rc, J, K, ncomp, nfront, nvtx, v ;
int      *bndwghts, *compids, *fch, *map, *nodwghts, 
         *par, *sib, *vtxToFront ;
IV       *compidsIV, *mapIV ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 5 ) {
   fprintf(stdout, 
      "\n\n usage : %s msglvl msgFile inETreeFile outIVfile"
      "\n    msglvl      -- message level"
      "\n    msgFile     -- message file"
      "\n    inETreeFile -- input file, must be *.etreef or *.etreeb"
      "\n    outIVfile   -- output file, must be *.ivf or *.ivb"
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
inETreeFileName = argv[3] ;
outIVfileName   = argv[4] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl      -- %d" 
        "\n msgFile     -- %s" 
        "\n inETreeFile -- %s" 
        "\n outIVfile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, outIVfileName) ;
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
fflush(msgFile) ;
nfront     = ETree_nfront(etree) ;
nvtx       = ETree_nvtx(etree) ;
bndwghts   = ETree_bndwghts(etree) ;
vtxToFront = ETree_vtxToFront(etree) ;
nodwghts   = ETree_nodwghts(etree) ;
par        = ETree_par(etree) ;
fch        = ETree_fch(etree) ;
sib        = ETree_sib(etree) ;
tree       = ETree_tree(etree) ;
/*
   -----------------------------------------
   create the map from fronts to components,
   top level separator(s) are component zero
   -----------------------------------------
*/
mapIV = IV_new() ;
IV_init(mapIV, nfront, NULL) ;
map = IV_entries(mapIV) ;
ncomp = 0 ;
for ( J = Tree_preOTfirst(tree) ;
      J != -1 ;
      J = Tree_preOTnext(tree, J) ) { 
   if ( (K = par[J]) == -1 ) {
      map[J] = 0 ;
   } else if ( map[K] != 0 ) {
      map[J] = map[K] ;
   } else if ( J == fch[K] && sib[J] == -1 
            && bndwghts[J] == nodwghts[K] + bndwghts[K] ) {
      map[J] = 0 ;
   } else {
      map[J] = ++ncomp ;
   }
}
fprintf(msgFile, "\n\n mapIV object") ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(mapIV, msgFile) ;
} else {
   IV_writeStats(mapIV, msgFile) ;
}
/*
   ----------------------------------------
   fill the map from vertices to components
   ----------------------------------------
*/
compidsIV = IV_new() ;
IV_init(compidsIV, nvtx, NULL) ;
compids = IV_entries(compidsIV) ;
for ( v = 0 ; v < nvtx ; v++ ) {
   compids[v] = map[vtxToFront[v]] ;
}
fprintf(msgFile, "\n\n compidsIV object") ;
if ( msglvl > 2 ) {
   IV_writeForHumanEye(compidsIV, msgFile) ;
} else {
   IV_writeStats(compidsIV, msgFile) ;
}
fflush(msgFile) ;
/*
   -----------------------
   write out the IV object
   -----------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(compidsIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outIVfileName) ;
}
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from IV_writeToFile(%p,%s)",
           rc, compidsIV, outIVfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
IV_free(mapIV) ;
IV_free(compidsIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
