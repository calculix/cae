/*  testStorage.c  */

#include "../../ETree.h"
#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ----------------------------------------------------
   read in an ETree object.
   read in a Graph object.
   get the symbolic factorization IVL object.
   compute the storage profiles for the general sparse,
   forward sparse and multifrontal methods.
   
   created -- 96oct03, cca
   ----------------------------------------------------
*/
{
char     *inETreeFileName, *inGraphFileName ;
double   elapsed, nfops1, t1, t2 ;
double   *FSvec, *GSvec, *MFvec, *backwardops, *forwardops, *vmetric ;
DV       *vmetricDV ;
Graph    *graph ;
int      J, msglvl, nfent1, nfind1, nfront, nleaves1, nnode1, rc ;
IVL      *symbfacIVL ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 5 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inGraphFile "
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inGraphFile  -- input file, must be *.graphf or *.graphb"
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
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inGraphFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inETreeFileName, inGraphFileName) ;
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
/*
   --------------------------------------------
   create the symbolic factorization IVL object
   --------------------------------------------
*/
symbfacIVL = SymbFac_initFromGraph(etree, graph) ;
fprintf(msgFile, 
        "\n\n symbolic factorization IVL object in old ordering") ;
if ( msglvl > 2 ) {
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
} else {
   IVL_writeStats(symbfacIVL, msgFile) ;
}
fflush(msgFile) ;
/*
   --------------------------
   get the operations profile
   --------------------------
*/
vmetricDV = ETree_backwardOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC,
                              graph->vwghts, symbfacIVL) ;
vmetric = DV_entries(vmetricDV) ;
backwardops = DVinit(nfront, 0.0) ;
elapsed = 0.0 ;
for ( J = Tree_postOTfirst(etree->tree) ;
      J != -1 ;
      J = Tree_postOTnext(etree->tree, J) ) {
   elapsed += vmetric[J] ;
   backwardops[J] = elapsed ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n sum of backward ops = %.0f",
           DV_sum(vmetricDV)) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n backward ops") ;
   DVfprintf(msgFile, nfront, backwardops) ;
}
DV_free(vmetricDV) ;
vmetricDV = ETree_forwardOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
vmetric   = DV_entries(vmetricDV) ;
forwardops = DVinit(nfront, 0.0) ;
elapsed = 0.0 ;
for ( J = Tree_postOTfirst(etree->tree) ;
      J != -1 ;
      J = Tree_postOTnext(etree->tree, J) ) {
   elapsed += vmetric[J] ;
   forwardops[J] = elapsed ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n\n sum of forward ops = %.0f",
           DV_sum(vmetricDV)) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n forward ops") ;
   DVfprintf(msgFile, nfront, forwardops) ;
}
DV_free(vmetricDV) ;
/*
   --------------------------------------
   get the general sparse storage profile
   --------------------------------------
*/
GSvec = DVinit(nfront, 0.0) ;
ETree_GSstorageProfile(etree, SPOOLES_SYMMETRIC,
                       symbfacIVL, graph->vwghts, GSvec) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n GSvec storage") ;
   DVfprintf(msgFile, nfront, GSvec) ;
}
/*
   --------------------------------------
   get the forward sparse storage profile
   --------------------------------------
*/
FSvec = DVinit(nfront, 0.0) ;
ETree_FSstorageProfile(etree, SPOOLES_SYMMETRIC, symbfacIVL, FSvec) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n FSvec storage") ;
   DVfprintf(msgFile, nfront, FSvec) ;
}
/*
   ------------------------------------
   get the multifrontal storage profile
   ------------------------------------
*/
MFvec = DVinit(nfront, 0.0) ;
ETree_MFstackProfile(etree, SPOOLES_SYMMETRIC, MFvec) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n MFvec storage") ;
   DVfprintf(msgFile, nfront, MFvec) ;
}
if ( msglvl > 0 ) {
   fprintf(msgFile, 
"\n %%   five columns of data"
"\n %%   backward-ops GS-storage forward-ops FS-storage MF-storage") ;
   fprintf(msgFile, "\n data = [ ...") ;
   for ( J = Tree_postOTfirst(tree) ;
         J != -1 ;
         J = Tree_postOTnext(tree, J) ) {
/*
      fprintf(msgFile, "\n %12.0f %12.0f %12.0f %12.0f %12.0f", 
          backwardops[J], FSvec[J], forwardops[J], MGSvec[J], Fvec[J]) ;
*/
      fprintf(msgFile, "\n %12.0f %12.4e %12.0f, %12.4e %12.4e", 
              backwardops[J], GSvec[J]/nfent1, 
              forwardops[J], FSvec[J]/nfent1, MFvec[J]/nfent1) ;
   }
   fprintf(msgFile, 
    " ] ;"
    "\n bops = data(:,1) ;"
    "\n gs   = data(:,2) ;"
    "\n fops = data(:,3) ;"
    "\n fs   = data(:,4) ;"
    "\n mf   = data(:,5) ;"
    "\n\n plot( bops, gs, '-o', fops, fs, '-v', fops, mf, '-s') ; "
    "\n xmax = max(bops) ;"
    "\n ymax = max( [ max(gs) max(fs) max(mf) ] ) ;"
    "\n axis([0, xmax, 0, ymax]) ;"
    "\n xlabel(' elapsed operations') ;"
    "\n ylabel(' fraction of total factor storage') ;"
    "\n title(' workspace profile, ""x"" GS, ""o"" FS, ""*"" MF') ;"
    "\n text( 0.1*xmax, 0.9*ymax, 'circle   -- general sparse') ;"
    "\n text( 0.1*xmax, 0.8*ymax, 'triangle -- forward sparse') ;"
    "\n text( 0.1*xmax, 0.7*ymax, 'square   -- multifrontal') ;" ) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
Graph_free(graph) ;
IVL_free(symbfacIVL) ;
DVfree(GSvec) ;
DVfree(MFvec) ;
DVfree(forwardops) ;
DVfree(backwardops) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(-1) ; }

/*--------------------------------------------------------------------*/
