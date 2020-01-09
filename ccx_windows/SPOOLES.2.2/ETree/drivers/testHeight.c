/*  testHeight.c  */

#include "../../ETree.h"
#include "../../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------------------
   read in an ETree object. compute the height of the tree 
   w.r.t. an out-of-core forward sparse factorization
   
   created -- 99jan07, cca
   -------------------------------------------------------
*/
{
char     *inETreeFileName ;
double   nfops1, t1, t2 ;
IV       *dmetricIV, *vmetricIV ;
int      maxdepth, maxnent, msglvl, nfent1, nfind1, 
         nfront, nleaves1, nnode1, rc ;
ETree    *etree ;
FILE     *msgFile ;
Tree     *tree ;

if ( argc != 4 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile "
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
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
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName) ;
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
   ----------------------
   get the height profile
   ----------------------
*/
dmetricIV = Tree_setDepthImetric(tree, vmetricIV) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n entries depth per front") ;
   IV_writeForHumanEye(dmetricIV, msgFile) ;
   fflush(msgFile) ;
}
maxdepth = IV_max(dmetricIV) ;
fprintf(msgFile, "\n\n max depth = %d, fraction of total = %8.3f", 
        maxdepth, ((double) maxdepth)/nfent1) ;
fflush(msgFile) ;
fprintf(msgFile, "\n\n STATS : %12d %12d %12d %8.3f %8.3f",
        nfent1, maxnent, maxdepth, ((double) maxdepth)/maxnent,
        ((double) maxdepth)/nfent1) ;
{
int   J ;
int   *depth    = IV_entries(dmetricIV) ;
int   *par      = Tree_par(tree) ;
int   *fch      = Tree_fch(tree) ;
int   *sib      = Tree_sib(tree) ;
int   *nodwghts = ETree_nodwghts(etree) ;
int   *bndwghts = ETree_bndwghts(etree) ;
fprintf(msgFile, "\n\n J  par  fch  sib  |J|  |bndJ|  depth/maxdepth") ;
for ( J = 0 ; J < nfront ; J++ ) {
   fprintf(msgFile, "\n %7d %7d %7d %7d %7d %7d %8.3f",
        J, par[J], fch[J], sib[J], nodwghts[J], bndwghts[J],
        ((double) depth[J])/maxdepth) ;
}
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
IV_free(vmetricIV) ;
IV_free(dmetricIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
