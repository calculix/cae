/*  testSymbfacInpMtx.c  */

#include "../SymbFac.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   ---------------------------------------------------------------
   (1) read in an ETree object.
   (2) permute the ETree object
   (3) read in a InpMtx object.
   (4) permute the InpMtx object
   (5) get the symbolic factorization IVL object.
   (6) optionally write the old-to-new IV object to a file
   
   created -- 96oct03, cca
   ---------------------------------------------------------------
*/
{
char     *inETreeFileName, *inInpMtxFileName, *outETreeFileName,
         *outIVfileName, *outIVLfileName ;
double   nfops1, t1, t2 ;
InpMtx    *inpmtx ;
int      msglvl, nfent1, nfind1, nleaves1, nnode1, rc ;
IV       *vtxOldToNewIV ;
IVL      *symbfacIVL ;
ETree    *etree ;
FILE     *msgFile ;

if ( argc != 8 ) {
   fprintf(stdout, 
"\n\n usage : %s msglvl msgFile inETreeFile inInpMtxFile outETreeFile"
"\n         outIVfile outIVLfile"
"\n    msglvl       -- message level"
"\n    msgFile      -- message file"
"\n    inETreeFile  -- input file, must be *.etreef or *.etreeb"
"\n    inInpMtxFile  -- input file, must be *.inpmtxf or *.inpmtxb"
"\n    outETreeFile -- output file, must be *.etreef or *.etreeb"
"\n    outIVfile    -- output file for oldToNew vector,"
"\n                    must be *.ivf or *.ivb"
"\n    outIVLfile   -- output file for symbolic factorization object"
"\n                    must be *.ivlf or *.ivlb"
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
inInpMtxFileName  = argv[4] ;
outETreeFileName = argv[5] ;
outIVfileName    = argv[6] ;
outIVLfileName   = argv[7] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl        -- %d" 
        "\n msgFile       -- %s" 
        "\n inETreeFile   -- %s" 
        "\n inInpMtxFile  -- %s" 
        "\n outETreeFile  -- %s" 
        "\n outIVfile     -- %s" 
        "\n outIVLfile    -- %s" 
        "\n",
        argv[0], msglvl, argv[2], 
        inETreeFileName, inInpMtxFileName, outETreeFileName, 
        outIVfileName, outIVLfileName) ;
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
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading ETree object from file %s",
           inETreeFileName) ;
   if ( msglvl == 2 ) {
      ETree_writeStats(etree, msgFile) ;
   } else {
      ETree_writeForHumanEye(etree, msgFile) ;
   }
}
fflush(msgFile) ;
/*
   ----------------------
   compute the statistics
   ----------------------
*/
MARKTIME(t1) ;
nnode1   = etree->tree->n ;
nfind1   = ETree_nFactorIndices(etree) ;
nfent1   = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
nfops1   = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
nleaves1 = Tree_nleaves(etree->tree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compute statistics", t2 - t1) ;
fprintf(msgFile, 
        "\n %d nodes, %d indices, %d entries, %.0f ops, %d leaves",
        nnode1, nfind1, nfent1, nfops1, nleaves1) ;
fflush(msgFile) ;
/*
   -----------------------------
   get the permutation IV object
   -----------------------------
*/
MARKTIME(t1) ;
vtxOldToNewIV = ETree_oldToNewVtxPerm(etree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : get permutation", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n vertex old-to-new IV object") ;
   if ( msglvl == 2 ) {
      IV_writeStats(vtxOldToNewIV, msgFile) ;
   } else {
      IV_writeForHumanEye(vtxOldToNewIV, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   ----------------------------------------
   permute the vertices in the ETree object
   ----------------------------------------
*/
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n before permuting the vertex map") ;
   if ( msglvl == 2 ) {
      ETree_writeStats(etree, msgFile) ;
   } else {
      ETree_writeForHumanEye(etree, msgFile) ;
   }
   fflush(msgFile) ;
}
MARKTIME(t1) ;
ETree_permuteVertices(etree, vtxOldToNewIV) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : permute ETree", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after permuting the vertex map") ;
   if ( msglvl == 2 ) {
      ETree_writeStats(etree, msgFile) ;
   } else {
      ETree_writeForHumanEye(etree, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   -------------------------
   read in the InpMtx object
   -------------------------
*/
if ( strcmp(inInpMtxFileName, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
inpmtx = InpMtx_new() ;
MARKTIME(t1) ;
rc = InpMtx_readFromFile(inpmtx, inInpMtxFileName) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : read in inpmtx from file %s",
        t2 - t1, inInpMtxFileName) ;
if ( rc != 1 ) {
   fprintf(msgFile, 
           "\n return value %d from InpMtx_readFromFile(%p,%s)",
           rc, inpmtx, inInpMtxFileName) ;
   exit(-1) ;
}
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after reading InpMtx object from file %s",
           inInpMtxFileName) ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(inpmtx, msgFile) ;
   } else {
      InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   }
   fflush(msgFile) ;
}
if ( INPMTX_IS_BY_ROWS(inpmtx) ) {
   fprintf(msgFile, "\n matrix coordinate type is rows") ;
} else if ( INPMTX_IS_BY_COLUMNS(inpmtx) ) {
   fprintf(msgFile, "\n matrix coordinate type is columns") ;
} else if ( INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   fprintf(msgFile, "\n matrix coordinate type is chevrons") ;
} else {
   fprintf(msgFile, "\n\n, error, bad coordinate type") ;
   exit(-1) ;
}
if ( INPMTX_IS_RAW_DATA(inpmtx) ) {
   fprintf(msgFile, "\n matrix storage mode is raw data\n") ;
} else if ( INPMTX_IS_SORTED(inpmtx) ) {
   fprintf(msgFile, "\n matrix storage mode is sorted\n") ;
} else if ( INPMTX_IS_BY_VECTORS(inpmtx) ) {
   fprintf(msgFile, "\n matrix storage mode is by vectors\n") ;
} else {
   fprintf(msgFile, "\n\n, error, bad storage mode") ;
   exit(-1) ;
}
/*
   --------------------------
   permute the InpMtx object  
   --------------------------
*/
MARKTIME(t1) ;
InpMtx_permute(inpmtx, IV_entries(vtxOldToNewIV), 
                       IV_entries(vtxOldToNewIV)) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : permute the matrix", t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n after permuting InpMtx object") ;
   if ( msglvl == 2 ) {
      InpMtx_writeStats(inpmtx, msgFile) ;
   } else {
      InpMtx_writeForHumanEye(inpmtx, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   ---------------------------------
   map entries to the upper triangle
   ---------------------------------
*/
MARKTIME(t1) ;
InpMtx_mapToUpperTriangle(inpmtx) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : map to upper triangle", t2 - t1) ;
/*
   -----------------------------------------------------
   change coordinate type and storage mode, if necessary
   -----------------------------------------------------
*/
if ( ! INPMTX_IS_BY_CHEVRONS(inpmtx) ) {
   MARKTIME(t1) ;
   InpMtx_changeCoordType(inpmtx, INPMTX_BY_CHEVRONS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : change coordinate type", t2 - t1) ;
}
if ( INPMTX_IS_RAW_DATA(inpmtx) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(inpmtx, INPMTX_SORTED) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : sort entries ", t2 - t1) ;
}
if ( INPMTX_IS_SORTED(inpmtx) ) {
   MARKTIME(t1) ;
   InpMtx_changeStorageMode(inpmtx, INPMTX_BY_VECTORS) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : convert to vectors ", t2 - t1) ;
}
/*
   --------------------------------------------
   create the symbolic factorization IVL object
   --------------------------------------------
*/
MARKTIME(t1) ;
symbfacIVL = SymbFac_initFromInpMtx(etree, inpmtx) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compute the symbolic factorization",
        t2 - t1) ;
if ( msglvl > 1 ) {
   fprintf(msgFile, "\n\n symbolic factorization IVL object") ;
   if ( msglvl == 2 ) {
      IVL_writeStats(symbfacIVL, msgFile) ;
   } else {
      IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   }
   fflush(msgFile) ;
}
/*
   --------------------------
   compute the new statistics
   --------------------------
*/
MARKTIME(t1) ;
nnode1   = etree->tree->n ;
nfind1   = ETree_nFactorIndices(etree) ;
nfent1   = ETree_nFactorEntries(etree, SPOOLES_SYMMETRIC) ;
nfops1   = ETree_nFactorOps(etree, SPOOLES_REAL, SPOOLES_SYMMETRIC) ;
nleaves1 = Tree_nleaves(etree->tree) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %9.5f : compute statistics", t2 - t1) ;
fprintf(msgFile, 
        "\n %d nodes, %d indices, %d entries, %.0f ops, %d leaves",
        nnode1, nfind1, nfent1, nfops1, nleaves1) ;
fflush(msgFile) ;
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
if ( rc != 1 ) {
   fprintf(msgFile, "\n return value %d from ETree_writeToFile(%p,%s)",
           rc, etree, outETreeFileName) ;
}
/*
   ----------------------------------------------
   optionally write out the permutation IV object
   ----------------------------------------------
*/
if ( strcmp(outIVfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IV_writeToFile(vtxOldToNewIV, outIVfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write permutation IV to file %s",
           t2 - t1, outIVfileName) ;
}
/*
   -----------------------------------
   optionally write out the IVL object
   -----------------------------------
*/
if ( strcmp(outIVLfileName, "none") != 0 ) {
   MARKTIME(t1) ;
   rc = IVL_writeToFile(symbfacIVL, outIVLfileName) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %9.5f : write symbfac IVL to file %s",
           t2 - t1, outIVLfileName) ;
}
/*
   ----------------
   free the objects
   ----------------
*/
ETree_free(etree) ;
InpMtx_free(inpmtx) ;
IV_free(vtxOldToNewIV) ;
IVL_free(symbfacIVL) ;

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
