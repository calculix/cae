/*  postponed.c  */

#include "../FrontMtx.h"

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------------
   purpose -- to assemble any postponed data into frontJ

   frontJ  -- pointer to Chv objec that contains current front
   chvlist -- pointer to a ChvList object that handles the
              lists of postponed Chv objects
   chvmanager -- pointer to a ChvManager object for the list
                 of free Chv objects
   pndelay -- pointer to address to contain the # of delayed rows 
              and columns that were assembled into the front

   return value -- pointer to Chv object that contains the new front

   created -- 98may04, cca
   ------------------------------------------------------------------
*/
Chv *
FrontMtx_assemblePostponedData (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   ChvList      *chvlist,
   ChvManager   *chvmanager,
   int          *pndelay
) {
Chv   *child, *child2, *firstchild, *newfrontJ, *nextchild, *prev ;
int    nbytes, nDnew ;

if ( (firstchild = ChvList_getList(chvlist, frontJ->id)) == NULL ) {
/*
   -------------------------------------
   quick return, no children to assemble
   -------------------------------------
*/
   *pndelay = 0 ;
   return(frontJ) ;
}
/*
   -------------------------------------------------------
   order the children in ascending order of their front id
   this is done to ensure that the serial, multithreaded
   and MPI codes all assemble the same frontal matrix.
   -------------------------------------------------------
*/
#if MYDEBUG > 0
fprintf(stdout, "\n postponed children of %d :", frontJ->id) ;
for ( child = firstchild ; child != NULL ; child = child->next ) {
   fprintf(stdout, " %d", child->id) ;
}
fflush(stdout) ;
#endif
for ( child = firstchild, firstchild = NULL ; 
      child != NULL ; 
      child = nextchild ) {
   nextchild = child->next ;
   for ( child2 = firstchild, prev = NULL ; 
         child2 != NULL ; 
         child2 = child2->next ) {
      if ( child2->id > child->id ) {
         break ;
      }
      prev = child2 ;
   }
   if ( prev == NULL ) {
      firstchild = child ;
   } else {
      prev->next = child ;
   }
   child->next = child2 ;
}
#if MYDEBUG > 0
fprintf(stdout, "\n front %d, postponed children reordered :", 
        frontJ->id) ;
for ( child = firstchild ; child != NULL ; child = child->next ) {
   fprintf(stdout, " %d", child->id) ;
}
fflush(stdout) ;
#endif
/*
   --------------------------
   compute the new dimensions
   --------------------------
*/
nDnew = frontJ->nD ;
for ( child = firstchild ; child != NULL ; child = child->next ) {
   nDnew += child->nD ;
}
/*
   ------------------------
   get a new chevron object
   ------------------------
*/
nbytes = Chv_nbytesNeeded(nDnew, frontJ->nL, frontJ->nU, 
                          frontJ->type, frontJ->symflag) ;
newfrontJ = ChvManager_newObjectOfSizeNbytes(chvmanager, nbytes) ;
Chv_init(newfrontJ, frontJ->id, nDnew, frontJ->nL, frontJ->nU, 
         frontJ->type, frontJ->symflag) ;
/*
   ----------------------------------------------------------
   pivoting has been enabled, assemble any postponed chevrons
   ----------------------------------------------------------
*/
*pndelay = Chv_assemblePostponedData(newfrontJ, frontJ, firstchild) ;
/*
   --------------------------------------------------
   now put the postponed chevrons onto the free list.
   --------------------------------------------------
*/
ChvManager_releaseListOfObjects(chvmanager, firstchild) ;
/*
   -------------------------------------
   set the delay to zero if a root front
   -------------------------------------
*/
if ( frontJ->nU == 0 ) {
   *pndelay = 0 ;
}
return(newfrontJ) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------
   purpose -- extract and store the postponed data

   frontJ  -- pointer to present front object
   npost   -- # of postponed rows and columns in frontJ
   K       -- parent of J
   chvlist -- pointer to a ChvList object that handles the
              lists of postponed Chv objects
              a singly linked list to assemble
   chvmanager -- pointer to a ChvManager object for the list
                 of free Chv objects

   created -- 98may04, cca
   ---------------------------------------------------------
*/
void
FrontMtx_storePostponedData (
   FrontMtx     *frontmtx,
   Chv          *frontJ,
   int          npost,
   int          K,
   ChvList      *chvlist,
   ChvManager   *chvmanager
) {
Chv   *chv ;
int    nbytes, nD, nent, nind, nL, nU ;

if ( npost <= 0 && chvlist != NULL ) {
   if ( K == -1 ) {
      ChvList_addObjectToList(chvlist, NULL, frontmtx->nfront) ;
   } else {
      ChvList_addObjectToList(chvlist, NULL, K) ;
   }
   return ;
} 
/*
   --------------------------------------
   find the number of indices and entries 
   necessary to store the delayed data
   --------------------------------------
*/
Chv_dimensions(frontJ, &nD, &nL, &nU) ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n front %d: npost = %d, nD = %d, nL = %d, nU = %d", 
        frontJ->id, npost, nD, nL, nU) ;
fflush(stdout) ;
#endif
if ( CHV_IS_SYMMETRIC(frontJ) || CHV_IS_HERMITIAN(frontJ) ) {
   nind = npost + nU ;
   nent = (npost*(npost+1))/2 + npost*nU ;
} else if ( CHV_IS_NONSYMMETRIC(frontJ) ) {
   nind = 2*(npost + nU) ;
   nent = npost*(npost + 2*nU) ;
}
/*
   ------------------------------------
   get a Chv object from the free list
   ------------------------------------
*/
nbytes = Chv_nbytesNeeded(npost, nL, nU, frontJ->type, frontJ->symflag);
chv = ChvManager_newObjectOfSizeNbytes(chvmanager, nbytes) ;
Chv_init(chv, frontJ->id, npost, nL, nU, frontJ->type, frontJ->symflag);
/*
   ----------------------
   store the delayed data
   ----------------------
*/
Chv_copyTrailingPortion(chv, frontJ, nD - npost) ;
frontJ->nD -= npost ;
frontJ->nL += npost ;
frontJ->nU += npost ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n postponed chevron %p", chv) ;
Chv_writeForHumanEye(chv, stdout) ;
fflush(stdout) ;
#endif
/*
   ------------------------------
   link the postponed Chv object
   ------------------------------
*/
if ( K == -1 ) {
   ChvList_addObjectToList(chvlist, chv, frontmtx->nfront) ;
} else {
   ChvList_addObjectToList(chvlist, chv, K) ;
}
return ; }
   
/*--------------------------------------------------------------------*/
