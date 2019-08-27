/*  Setup.c  */

#include "../Bridge.h"

#define MYDEBUG 1

#if MYDEBUG > 0
static int count_Setup = 0 ;
static double time_Setup = 0.0 ;
#endif

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------------
   purpose --

   given InpMtx objects that contain A and B, initialize the bridge
   data structure for the serial factor's, solve's and mvm's.

   data -- pointer to a Bridge object
   pprbtype -- pointer to value containing problem type
     *prbtype = 1 --> A X = B X Lambda, vibration problem
     *prbtype = 2 --> A X = B X Lambda, buckling problem
     *prbtype = 3 --> A X = X Lambda, simple eigenvalue problem
   pneqns  -- pointer to value containing number of equations
   pmxbsz  -- pointer to value containing blocksize
   A       -- pointer to InpMtx object containing A
   B       -- pointer to InpMtx object containing B
   pseed   -- pointer to value containing a random number seed
   pmsglvl -- pointer to value containing a message level
   msgFile -- message file pointer

   return value --
      1 -- normal return
     -1 -- data is NULL
     -2 -- pprbtype is NULL
     -3 -- *pprbtype is invalid
     -4 -- pneqns is NULL
     -5 -- *pneqns is invalid
     -6 -- pmxbsz is NULL
     -7 -- *pmxbsz is invalid
     -8 -- A and B are NULL
     -9 -- pseed is NULL
    -10 -- pmsglvl is NULL
    -11 -- *pmsglvl > 0 and msgFile is NULL

   created -- 98aug10, cca
   ----------------------------------------------------------------
*/
int
Setup (
   void     *data,
   int      *pprbtype,
   int      *pneqns,
   int      *pmxbsz,
   InpMtx   *A,
   InpMtx   *B,
   int      *pseed,
   int      *pmsglvl,
   FILE     *msgFile
) {
Bridge   *bridge = (Bridge *) data ;
double   sigma[2] ;
Graph    *graph ;
int      maxdomainsize, maxsize, maxzeros, msglvl, mxbsz, 
         nedges, neqns, prbtype, seed ;
IVL      *adjIVL ;
#if MYDEBUG > 0
double   t1, t2 ;
MARKTIME(t1) ;
count_Setup++ ;
fprintf(stdout, "\n (%d) Setup()", count_Setup) ;
fflush(stdout) ;
#endif
/*
   --------------------
   check the input data
   --------------------
*/
if ( data == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n data is NULL\n") ;
   return(-1) ;
}
if ( pprbtype == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n prbtype is NULL\n") ;
   return(-2) ;
}
prbtype = *pprbtype ;
if ( prbtype < 1 || prbtype > 3 ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n prbtype = %d, is invalid\n", prbtype) ;
   return(-3) ;
}
if ( pneqns == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n pneqns is NULL\n") ;
   return(-4) ;
}
neqns = *pneqns ;
if ( neqns <= 0 ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n neqns = %d, is invalid\n", neqns) ;
   return(-5) ;
}
if ( pmxbsz == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n pmxbsz is NULL\n") ;
   return(-6) ;
}
mxbsz = *pmxbsz ;
if ( mxbsz <= 0 ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n *pmxbsz = %d, is invalid\n", mxbsz) ;
   return(-7) ;
}
if ( A == NULL && B == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n A and B are NULL\n") ;
   return(-8) ;
}
if ( pseed == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n pseed is NULL\n") ;
   return(-9) ;
}
seed = *pseed ;
if ( pmsglvl == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n pmsglvl is NULL\n") ;
   return(-10) ;
}
msglvl = *pmsglvl ;
if ( msglvl > 0 && msgFile == NULL ) {
   fprintf(stderr, "\n fatal error in Setup()"
           "\n msglvl = %d, msgFile = NULL\n", msglvl) ;
   return(-11) ;
}
bridge->msglvl  = msglvl  ;
bridge->msgFile = msgFile ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n inside Setup()"
           "\n neqns = %d, prbtype = %d, mxbsz = %d, seed = %d",
           neqns, prbtype, mxbsz, seed) ;
   if ( A != NULL ) {
      fprintf(msgFile, "\n\n matrix A") ;
      InpMtx_writeForHumanEye(A, msgFile) ;
   }
   if ( B != NULL ) {
      fprintf(msgFile, "\n\n matrix B") ;
      InpMtx_writeForHumanEye(B, msgFile) ;
   }
   fflush(msgFile) ;
}
bridge->prbtype = prbtype ;
bridge->neqns   = neqns   ;
bridge->mxbsz   = mxbsz   ;
bridge->A       = A       ;
bridge->B       = B       ;
bridge->seed    = seed    ;
/*
   ----------------------------
   create and initialize pencil
   ----------------------------
*/
sigma[0] = 1.0; sigma[1] = 0.0;
bridge->pencil = Pencil_new() ;
Pencil_setDefaultFields(bridge->pencil) ;
Pencil_init(bridge->pencil, SPOOLES_REAL, SPOOLES_SYMMETRIC,
            A, sigma, B) ;
/*
   --------------------------------
   convert to row or column vectors
   --------------------------------
*/
if ( A != NULL ) {
   if ( ! INPMTX_IS_BY_ROWS(A) && ! INPMTX_IS_BY_COLUMNS(A) ) {
      InpMtx_changeCoordType(A, INPMTX_BY_ROWS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(A) ) {
      InpMtx_changeStorageMode(A, INPMTX_BY_VECTORS) ;
   }
}
if ( B != NULL ) {
   if ( ! INPMTX_IS_BY_ROWS(B) && ! INPMTX_IS_BY_COLUMNS(B) ) {
      InpMtx_changeCoordType(B, INPMTX_BY_ROWS) ;
   }
   if ( ! INPMTX_IS_BY_VECTORS(B) ) {
      InpMtx_changeStorageMode(B, INPMTX_BY_VECTORS) ;
   }
}
/*
   -------------------------------
   create a Graph object for A + B
   -------------------------------
*/
graph  = Graph_new() ;
adjIVL = Pencil_fullAdjacency(bridge->pencil) ;
nedges = IVL_tsize(adjIVL),
Graph_init2(graph, 0, bridge->neqns, 0, nedges,
            bridge->neqns, nedges, adjIVL, NULL, NULL) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n graph of the input matrix") ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
/*
   ---------------
   order the graph
   ---------------
*/
maxdomainsize = neqns / 64 ;
if ( maxdomainsize == 0 ) {
   maxdomainsize = 1 ;
}
maxzeros  = (int) (0.01*neqns) ;
maxsize   = 64 ;
bridge->frontETree = orderViaBestOfNDandMS(graph, maxdomainsize, 
                     maxzeros, maxsize, bridge->seed, msglvl, msgFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree from ordering") ;
   ETree_writeForHumanEye(bridge->frontETree, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------------------
   get the old-to-new and new-to-old permutations
   ----------------------------------------------
*/
bridge->oldToNewIV = ETree_oldToNewVtxPerm(bridge->frontETree) ;
bridge->newToOldIV = ETree_newToOldVtxPerm(bridge->frontETree) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n old-to-new permutation") ;
   IV_writeForHumanEye(bridge->oldToNewIV, msgFile) ;
   fprintf(msgFile, "\n\n new-to-old permutation") ;
   IV_writeForHumanEye(bridge->newToOldIV, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------
   permute the vertices in the front tree
   --------------------------------------
*/
ETree_permuteVertices(bridge->frontETree, bridge->oldToNewIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted front etree") ;
   ETree_writeForHumanEye(bridge->frontETree, msgFile) ;
   fflush(msgFile) ;
}
/*
   -------------------------------------------
   permute the entries in the pencil. 
   note, after the permutation the 
   entries are mapped into the upper triangle.
   -------------------------------------------
*/
Pencil_permute(bridge->pencil, bridge->oldToNewIV, bridge->oldToNewIV) ;
Pencil_mapToUpperTriangle(bridge->pencil) ;
Pencil_changeCoordType(bridge->pencil, INPMTX_BY_CHEVRONS) ;
Pencil_changeStorageMode(bridge->pencil, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted pencil") ;
   Pencil_writeForHumanEye(bridge->pencil, msgFile) ;
   fflush(msgFile) ;
}
/*
   ----------------------------------
   compute the symbolic factorization
   ----------------------------------
*/
bridge->symbfacIVL = SymbFac_initFromPencil(bridge->frontETree, 
                                            bridge->pencil) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n symbolic factorization") ;
   IVL_writeForHumanEye(bridge->symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*
   --------------------------------------------------
   create a FrontMtx object to hold the factorization
   --------------------------------------------------
*/
bridge->frontmtx = FrontMtx_new() ;
/*
   ------------------------------------------------------------
   create a SubMtxManager object to hold the factor submatrices
   ------------------------------------------------------------
*/
bridge->mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(bridge->mtxmanager, NO_LOCK, 0) ;
/*
   ------------------------------------------------------------
   allocate the working objects X and Y for the matrix multiply
   ------------------------------------------------------------
*/
bridge->X = DenseMtx_new() ;
DenseMtx_init(bridge->X, SPOOLES_REAL, 0, 0, neqns, mxbsz, 1, neqns) ;
bridge->Y = DenseMtx_new() ;
DenseMtx_init(bridge->Y, SPOOLES_REAL, 0, 0, neqns, mxbsz, 1, neqns) ;
/*
   ------------------------
   free the working storage
   ------------------------
*/
Graph_free(graph) ;

#if MYDEBUG > 0
MARKTIME(t2) ;
time_Setup += t2 - t1 ;
fprintf(stdout, ", %8.3f seconds, %8.3f total time", 
        t2 - t1, time_Setup) ;
fflush(stdout) ;
#endif

return(1) ; }

/*--------------------------------------------------------------------*/
