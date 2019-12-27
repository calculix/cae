/*  permuteETree.c  */

#include "../ETree.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -----------------------------------------------------------
   read in an ETree object.
   if the ETree is defined on a compressed graph
      read in an equivalence map IV object.
      expand the ETree to be defined on the unit weight graph.
   endif
   get the old-to-new vertex permutation.
   permute the vtx-to-front map.

   created -- 97feb28, cca
   -----------------------------------------------------------
*/
{
char     *inEqmapIVfileName, *inETreeFileName, *outETreeFileName,
         *outIVfileName ;
double   t1, t2 ;
int      msglvl, rc ;
IV       *eqmapIV, *vtxOldToNewIV ;
ETree    *etree, *etree2 ;
FILE     *msgFile ;

if ( argc != 7 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inEqmapFile "
"\n         outETreeFile outIVfile "
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inEqmapFile  -- input file, must be *.ivf or *.ivb"
"\n    outETreeFile -- output file, must be *.etreef or *.etreeb"
"\n    outIVfile    -- output file for oldToNew vector,"
"\n                    must be *.ivf or *.ivb"
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
inEqmapIVfileName = argv[4] ;
outETreeFileName  = argv[5] ;
outIVfileName     = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inEqmapFile   -- %s" 
        "\n outETreeFile  -- %s" 
        "\n outIVfile     -- %s" 
        "\n",
        argv[0], msglvl, argv[2], inETreeFileName, inEqmapIVfileName, 
        outETreeFileName, outIVfileName) ;
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

if ( strcmp(inEqmapIVfileName, "none") != 0 ) {
/*
   -------------------------------------
   read in the equivalence map IV object
   -------------------------------------
*/
   eqmapIV = IV_new() ;
   MARKTIME(t1) ;
   rc = IV_readFromFile(eqmapIV, inEqmapIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : read in eqmapIV from file %s",
           t2 - t1, inEqmapIVfileName) ;
   if ( rc != 1 ) {
      fprintf(msgFile, "\n return value %d from IV_readFromFile(%p,%s)",
              rc, eqmapIV, inEqmapIVfileName) ;
      exit(-1) ;
   }
   fprintf(msgFile, "\n\n after reading IV object from file %s",
           inEqmapIVfileName) ;
   if ( msglvl > 2 ) {
      IV_writeForHumanEye(eqmapIV, msgFile) ;
   } else {
      IV_writeStats(eqmapIV, msgFile) ;
   }
   fflush(msgFile) ;
/*
   ---------------------
   expand the front tree
   ---------------------
*/
   MARKTIME(t1) ;
   etree2 = ETree_expand(etree, eqmapIV) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n\n CPU %9.5f : expand the ETree", t2 - t1) ;
   fprintf(msgFile, "\n\n expanded ETree") ;
   if ( msglvl > 2 ) {
      ETree_writeForHumanEye(etree2, msgFile) ;
   } else {
      ETree_writeStats(etree2, msgFile) ;
   }
/*
   ------------------------------------------------
   free the old ETree object and the eqmapIV object
   ------------------------------------------------
*/
   ETree_free(etree) ;
   etree = etree2 ;
   etree2 = NULL ;
   IV_free(eqmapIV) ;
}
/*
   -----------------------------
   get the permutation IV object
   -----------------------------
*/
MARKTIME(t1) ;
vtxOldToNewIV = ETree_oldToNewVtxPerm(etree) ;
MARKTIME(t2) ;
fprintf(msgFile, 
        "\n\n CPU %9.5f : get the old-to-new permutation", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   IV_writeForHumanEye(vtxOldToNewIV, msgFile) ;
} else {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   IV_writeStats(vtxOldToNewIV, msgFile) ;
}
fflush(msgFile) ;
/*
   ------------------------------------------------
   overwrite the ETree object with the new ordering
   ------------------------------------------------
*/
ETree_permuteVertices(etree, vtxOldToNewIV) ;
fprintf(msgFile, "\n\n after permuting the vertices") ;
if ( msglvl > 2 ) {
   ETree_writeForHumanEye(etree, msgFile) ;
} else {
   ETree_writeStats(etree, msgFile) ;
}
/*
   -------------------------------------
   optionally write out the ETree object
   -------------------------------------
*/
if ( strcmp(outETreeFileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = ETree_writeToFile(etree, outETreeFileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write etree to file %s",
           t2 - t1, outETreeFileName) ;
}
/*
   ----------------------------------
   optionally write out the IV object
   ----------------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(vtxOldToNewIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write vtxOldToNewIV to file %s",
           t2 - t1, outIVfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
IV_free(vtxOldToNewIV) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
