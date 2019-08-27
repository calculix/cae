/*  IO.c  */

#include "../FrontMtx.h"

static const char *suffixb = ".frontmtxb" ;
static const char *suffixf = ".frontmtxf" ;

#define MYDEBUG 0

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   purpose -- to read in an FrontMtx object from a file

   input --

      fn -- filename, must be *.frontmtxb or *.frontmtxf

   return value -- 1 if success, 0 if failure

   created -- 98may04, cca
   -----------------------------------------------------
*/
int
FrontMtx_readFromFile ( 
   FrontMtx   *frontmtx, 
   char       *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_readFromFile(%p,%s)"
           "\n bad input\n", frontmtx, fn) ;
   return(0) ;
}
/*
   -------------
   read the file
   -------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "rb")) == NULL ) {
         fprintf(stderr, "\n error in FrontMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", frontmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = FrontMtx_readFromBinaryFile(frontmtx, fp) ;
#if MYDEBUG > 0
         fprintf(stdout, 
                 "\n rc = %d from FrontMtx_readFromBinaryFile()", rc) ;
         fflush(stdout) ; 
#endif
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in FrontMtx_readFromFile(%p,%s)"
                 "\n unable to open file %s", frontmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = FrontMtx_readFromFormattedFile(frontmtx, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in FrontMtx_readFromFile(%p,%s)"
              "\n bad FrontMtx file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              frontmtx, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in FrontMtx_readFromFile(%p,%s)"
       "\n bad FrontMtx file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       frontmtx, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
#if MYDEBUG > 0
fprintf(stdout, 
        "\n returning rc = %d from FrontMtx_readFromFile()", rc) ;
fflush(stdout) ; 
#endif
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   purpose -- to read an FrontMtx object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 98may04, cca
   ------------------------------------------------------------
*/
int
FrontMtx_readFromFormattedFile ( 
   FrontMtx   *frontmtx, 
   FILE        *fp 
) {
SubMtx      *mtx ;
int       imtx, J, JK, KJ, nfront, nmtx, rc ;
int       itemp[10] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
FrontMtx_clearData(frontmtx) ;
/*
   -----------------------------
   read in the ten scalar fields
   -----------------------------
*/
if ( (rc = IVfscanf(fp, 10, itemp)) != 10 ) {
   fprintf(stderr, "\n error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", frontmtx, fp, rc, 10) ;
   return(0) ;
}
frontmtx->nfront       = nfront = itemp[0] ;
frontmtx->neqns        = itemp[1] ;
frontmtx->type         = itemp[2] ;
frontmtx->symmetryflag = itemp[3] ;
frontmtx->pivotingflag = itemp[4] ;
frontmtx->sparsityflag = itemp[5] ;
frontmtx->dataMode     = itemp[6] ;
frontmtx->nentD        = itemp[7] ;
frontmtx->nentL        = itemp[8] ;
frontmtx->nentU        = itemp[9] ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n\n nfront        = %d"
        "\n neqns         = %d"
        "\n type          = %d"
        "\n symmetryflag  = %d"
        "\n pivotingflag  = %d"
        "\n sparsityflag  = %d"
        "\n dataMode      = %d"
        "\n nentD         = %d"
        "\n nentL         = %d"
        "\n nentU         = %d",
        frontmtx->nfront, frontmtx->neqns, frontmtx->type,
        frontmtx->symmetryflag, frontmtx->pivotingflag, 
        frontmtx->sparsityflag, frontmtx->dataMode,
        frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
fflush(stdout) ;
#endif
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
frontmtx->frontETree = ETree_new() ;
rc = ETree_readFromFormattedFile(frontmtx->frontETree, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n error %d reading in frontETree object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
frontmtx->tree = frontmtx->frontETree->tree ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n ETree object") ;
ETree_writeForHumanEye(frontmtx->frontETree, stdout) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------
   read in the symbolic factorization object
   -----------------------------------------
*/
frontmtx->symbfacIVL = IVL_new() ;
frontmtx->symbfacIVL->type = IVL_CHUNKED ;
rc = IVL_readFromFormattedFile(frontmtx->symbfacIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
        "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
        "\n error %d reading in symbfacIVL object\n",
        frontmtx, fp, rc) ;
   return(0) ;
}
frontmtx->frontsizesIV = IV_new() ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   ---------------------------------
   read in the front sizes IV object
   ---------------------------------
*/
   rc = IV_readFromFormattedFile(frontmtx->frontsizesIV, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
           "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n error %d reading in frontsizesIV object\n",
           frontmtx, fp, rc) ;
      return(0) ;
   }
} else {
/*
   ----------------------------
   create from the ETree object
   ----------------------------
*/
   IV_init(frontmtx->frontsizesIV, frontmtx->nfront, NULL) ;
   IVcopy(nfront, IV_entries(frontmtx->frontsizesIV),
          ETree_nodwghts(frontmtx->frontETree)) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   -----------------------------------------
   read in the rowids and colids IVL objects
   -----------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->rowadjIVL = IVL_new() ;
      frontmtx->rowadjIVL->type = IVL_CHUNKED ;
      rc = IVL_readFromFormattedFile(frontmtx->rowadjIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in rowadjIVL object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   frontmtx->coladjIVL = IVL_new() ;
   frontmtx->coladjIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromFormattedFile(frontmtx->coladjIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in coladjIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
   fprintf(stdout, "\n coladjIVL read in") ;
   fflush(stdout) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
/*
   -------------------------------
   set up the five pointer vectors
   -------------------------------
*/
   ALLOCATE(frontmtx->p_mtxDJJ, struct _SubMtx *, nfront) ;
   ALLOCATE(frontmtx->p_mtxUJJ, struct _SubMtx *, nfront) ;
   ALLOCATE(frontmtx->p_mtxUJN, struct _SubMtx *, nfront) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      frontmtx->p_mtxDJJ[J] = NULL ;
      frontmtx->p_mtxUJJ[J] = NULL ;
      frontmtx->p_mtxUJN[J] = NULL ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      ALLOCATE(frontmtx->p_mtxLJJ, struct _SubMtx *, nfront) ;
      ALLOCATE(frontmtx->p_mtxLNJ, struct _SubMtx *, nfront) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         frontmtx->p_mtxLJJ[J] = NULL ;
         frontmtx->p_mtxLNJ[J] = NULL ;
      }
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -----------------------------
      read in the lower matrices
      -----------------------------
*/
      fscanf(fp, " %d", &nmtx) ;
      for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
         mtx = SubMtx_new() ;
         rc  = SubMtx_readFromFormattedFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
            return(0) ;
         }
         frontmtx->p_mtxLJJ[mtx->colid] = mtx ;
      }
      fscanf(fp, " %d", &nmtx) ;
      for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
         mtx = SubMtx_new() ;
         rc  = SubMtx_readFromFormattedFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
            return(0) ;
         }
         frontmtx->p_mtxLNJ[mtx->colid] = mtx ;
      }
   }
/*
   -----------------------------
   read in the diagonal matrices
   -----------------------------
*/
   fscanf(fp, " %d", &nmtx) ;
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromFormattedFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxDJJ[mtx->rowid] = mtx ;
   }
/*
   --------------------------
   read in the upper matrices
   --------------------------
*/
   fscanf(fp, " %d", &nmtx) ;
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc  = SubMtx_readFromFormattedFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
           "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n error %d reading in diag SubMtx object\n",
           frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxUJJ[mtx->rowid] = mtx ;
   }
   fscanf(fp, " %d", &nmtx) ;
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc  = SubMtx_readFromFormattedFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
           "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n error %d reading in diag SubMtx object\n",
           frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxUJN[mtx->rowid] = mtx ;
   }
} else {
/*
   -------------------------------------------------------
   read in the lower and upper block adjacency IVL objects
   -------------------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->lowerblockIVL = IVL_new() ;
      frontmtx->lowerblockIVL->type = IVL_CHUNKED ;
      rc = IVL_readFromFormattedFile(frontmtx->lowerblockIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in lowerblockIVL object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   frontmtx->upperblockIVL = IVL_new() ;
   frontmtx->upperblockIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromFormattedFile(frontmtx->upperblockIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in upperblockIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      --------------------------------------------------------------
      read in the lower submatrices and store them into a hash table
      --------------------------------------------------------------
*/
      frontmtx->lowerhash = I2Ohash_new() ;
      fscanf(fp, " %d", &nmtx) ;
      I2Ohash_init(frontmtx->lowerhash, nfront - 1, nmtx, 0) ;
      for ( KJ = 0 ; KJ < nmtx ; KJ++ ) {
         mtx = SubMtx_new() ;
         rc = SubMtx_readFromFormattedFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in lower SubMtx object\n",
              frontmtx, fp, rc) ;
            return(0) ;
         }
         I2Ohash_insert(frontmtx->lowerhash, 
                        mtx->rowid, mtx->colid, (void *) mtx) ;
      }
   }
/*
   -----------------------------------------------------------------
   read in the diagonal submatrices and store them into a hash table
   -----------------------------------------------------------------
*/
   ALLOCATE(frontmtx->p_mtxDJJ, struct _SubMtx *, nfront) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      frontmtx->p_mtxDJJ[J] = NULL ;
   }
   fscanf(fp, " %d", &nmtx) ;
   for ( J = 0 ; J < nmtx ; J++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromFormattedFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
              "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxDJJ[mtx->rowid] = mtx ;
   }
/*
   --------------------------------------------------------------
   read in the upper submatrices and store them into a hash table
   --------------------------------------------------------------
*/
   frontmtx->upperhash = I2Ohash_new() ;
   fscanf(fp, " %d", &nmtx) ;
   I2Ohash_init(frontmtx->upperhash, nfront-1, nmtx, 0) ;
   for ( JK = 0 ; JK < nmtx ; JK++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromFormattedFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
           "\n fatal error in FrontMtx_readFromFormattedFile(%p,%p)"
           "\n error %d reading in upper SubMtx object\n",
           frontmtx, fp, rc) ;
         return(0) ;
      }
      I2Ohash_insert(frontmtx->upperhash, 
                     mtx->rowid, mtx->colid, (void *) mtx) ;
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read an FrontMtx object from a binary file

   return value -- 1 if success, 0 if failure

   created -- 98may04, cca
   --------------------------------------------------------
*/
int
FrontMtx_readFromBinaryFile ( 
   FrontMtx   *frontmtx, 
   FILE        *fp 
) {
SubMtx      *mtx ;
int       imtx, J, JK, KJ, nfront, nmtx, rc ;
int       itemp[10] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
FrontMtx_clearData(frontmtx) ;
/*
   -----------------------------
   read in the ten scalar fields
   -----------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 10, fp)) != 10 ) {
   fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n %d items of %d read\n", frontmtx, fp, rc, 10) ;
   return(0) ;
}
frontmtx->nfront       = nfront = itemp[0] ;
frontmtx->neqns        = itemp[1] ;
frontmtx->type         = itemp[2] ;
frontmtx->symmetryflag = itemp[3] ;
frontmtx->pivotingflag = itemp[4] ;
frontmtx->sparsityflag = itemp[5] ;
frontmtx->dataMode     = itemp[6] ;
frontmtx->nentD        = itemp[7] ;
frontmtx->nentL        = itemp[8] ;
frontmtx->nentU        = itemp[9] ;
#if MYDEBUG > 0
fprintf(stdout, 
        "\n\n nfront        = %d"
        "\n neqns         = %d"
        "\n type          = %d"
        "\n symmetryflag  = %d"
        "\n pivotingflag  = %d"
        "\n sparsityflag  = %d"
        "\n dataMode      = %d"
        "\n nentD         = %d"
        "\n nentL         = %d"
        "\n nentU         = %d",
        frontmtx->nfront, frontmtx->neqns, frontmtx->type,
        frontmtx->symmetryflag, frontmtx->pivotingflag, 
        frontmtx->sparsityflag, frontmtx->dataMode,
        frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
fflush(stdout) ;
#endif
/*
   ------------------------
   read in the ETree object
   ------------------------
*/
frontmtx->frontETree = ETree_new() ;
rc = ETree_readFromBinaryFile(frontmtx->frontETree, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n error %d reading in frontETree object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
frontmtx->tree = frontmtx->frontETree->tree ;
#if MYDEBUG > 0
fprintf(stdout, "\n\n ETree object") ;
ETree_writeForHumanEye(frontmtx->frontETree, stdout) ;
fflush(stdout) ;
#endif
/*
   -----------------------------------------
   read in the symbolic factorization object
   -----------------------------------------
*/
frontmtx->symbfacIVL = IVL_new() ;
frontmtx->symbfacIVL->type = IVL_CHUNKED ;
rc = IVL_readFromBinaryFile(frontmtx->symbfacIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
        "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
        "\n error %d reading in symbfacIVL object\n",
        frontmtx, fp, rc) ;
   return(0) ;
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n symbfacIVL object") ;
   IVL_writeForHumanEye(frontmtx->symbfacIVL, stdout) ;
   fflush(stdout) ;
#endif

frontmtx->frontsizesIV = IV_new() ;
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   ---------------------------------
   read in the front sizes IV object
   ---------------------------------
*/
   rc = IV_readFromBinaryFile(frontmtx->frontsizesIV, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
           "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n error %d reading in frontsizesIV object\n",
           frontmtx, fp, rc) ;
      return(0) ;
   }
} else {
/*
   ----------------------------
   create from the ETree object
   ----------------------------
*/
   IV_init(frontmtx->frontsizesIV, frontmtx->nfront, NULL) ;
   IVcopy(nfront, IV_entries(frontmtx->frontsizesIV),
          ETree_nodwghts(frontmtx->frontETree)) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   -----------------------------------------
   read in the rowids and colids IVL objects
   -----------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->rowadjIVL = IVL_new() ;
      frontmtx->rowadjIVL->type = IVL_CHUNKED ;
      rc = IVL_readFromBinaryFile(frontmtx->rowadjIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
                 "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n error %d reading in rowadjIVL object\n",
                 frontmtx, fp, rc) ;
         return(0) ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n rowadjIVL object") ;
      IVL_writeForHumanEye(frontmtx->rowadjIVL, stdout) ;
      fflush(stdout) ;
#endif
   }
   frontmtx->coladjIVL = IVL_new() ;
   frontmtx->coladjIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromBinaryFile(frontmtx->coladjIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in coladjIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n coladjIVL object") ;
   IVL_writeForHumanEye(frontmtx->coladjIVL, stdout) ;
   fflush(stdout) ;
#endif
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
/*
   -------------------------------
   set up the five pointer vectors
   -------------------------------
*/
   ALLOCATE(frontmtx->p_mtxDJJ, struct _SubMtx *, nfront) ;
   ALLOCATE(frontmtx->p_mtxUJJ, struct _SubMtx *, nfront) ;
   ALLOCATE(frontmtx->p_mtxUJN, struct _SubMtx *, nfront) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      frontmtx->p_mtxDJJ[J] = NULL ;
      frontmtx->p_mtxUJJ[J] = NULL ;
      frontmtx->p_mtxUJN[J] = NULL ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      ALLOCATE(frontmtx->p_mtxLJJ, struct _SubMtx *, nfront) ;
      ALLOCATE(frontmtx->p_mtxLNJ, struct _SubMtx *, nfront) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         frontmtx->p_mtxLJJ[J] = NULL ;
         frontmtx->p_mtxLNJ[J] = NULL ;
      }
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -----------------------------
      read in the lower matrices
      -----------------------------
*/
      if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
         fprintf(stderr, 
                 "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
         return(0) ;
      }
      for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
         mtx = SubMtx_new() ;
         rc  = SubMtx_readFromBinaryFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
            return(0) ;
         }
         frontmtx->p_mtxLJJ[mtx->colid] = mtx ;
      }
      if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
         fprintf(stderr, 
                 "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
         return(0) ;
      }
      for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
         mtx = SubMtx_new() ;
         rc  = SubMtx_readFromBinaryFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
            return(0) ;
         }
         frontmtx->p_mtxLNJ[mtx->colid] = mtx ;
      }
   }
/*
   -----------------------------
   read in the diagonal matrices
   -----------------------------
*/
   if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
      return(0) ;
   }
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromBinaryFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxDJJ[mtx->rowid] = mtx ;
   }
/*
   --------------------------
   read in the upper matrices
   --------------------------
*/
   if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
      return(0) ;
   }
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc  = SubMtx_readFromBinaryFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
           "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n error %d reading in diag SubMtx object\n",
           frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxUJJ[mtx->rowid] = mtx ;
   }
   if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
      return(0) ;
   }
   for ( imtx = 0 ; imtx < nmtx ; imtx++ ) {
      mtx = SubMtx_new() ;
      rc  = SubMtx_readFromBinaryFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
           "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
           "\n error %d reading in diag SubMtx object\n",
           frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxUJN[mtx->rowid] = mtx ;
   }
} else {
/*
   -------------------------------------------------------
   read in the lower and upper block adjacency IVL objects
   -------------------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      frontmtx->lowerblockIVL = IVL_new() ;
      frontmtx->lowerblockIVL->type = IVL_CHUNKED ;
      rc = IVL_readFromBinaryFile(frontmtx->lowerblockIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
                 "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n error %d reading in lowerblockIVL object\n",
                 frontmtx, fp, rc) ;
         return(0) ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n lowerblockIVL object") ;
      IVL_writeForHumanEye(frontmtx->lowerblockIVL, stdout) ;
      fflush(stdout) ;
#endif
   }
   frontmtx->upperblockIVL = IVL_new() ;
   frontmtx->upperblockIVL->type = IVL_CHUNKED ;
   rc = IVL_readFromBinaryFile(frontmtx->upperblockIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in upperblockIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n upperblockIVL object") ;
   IVL_writeForHumanEye(frontmtx->upperblockIVL, stdout) ;
   fflush(stdout) ;
#endif
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      --------------------------------------------------------------
      read in the lower submatrices and store them into a hash table
      --------------------------------------------------------------
*/
      frontmtx->lowerhash = I2Ohash_new() ;
      if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
         fprintf(stderr, 
                 "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
         return(0) ;
      }
#if MYDEBUG > 0
      fprintf(stdout, "\n\n nmtx = %d", nmtx) ;
      fflush(stdout) ;
#endif
      I2Ohash_init(frontmtx->lowerhash, frontmtx->nfront-1, nmtx, 0) ;
      for ( KJ = 0 ; KJ < nmtx ; KJ++ ) {
         mtx = SubMtx_new() ;
         rc = SubMtx_readFromBinaryFile(mtx, fp) ;
         if ( rc != 1 ) {
            fprintf(stderr,
                 "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
                 "\n error %d reading in lower SubMtx object\n",
                 frontmtx, fp, rc) ;
            return(0) ;
         }
         I2Ohash_insert(frontmtx->lowerhash, 
                        mtx->rowid, mtx->colid, (void *) mtx) ;
#if MYDEBUG > 0
         fprintf(stdout, "\n\n --- lower SubMtx object") ;
         SubMtx_writeForHumanEye(mtx, stdout) ;
         fflush(stdout) ;
#endif
      }
   }
/*
   -----------------------------------------------------------------
   read in the diagonal submatrices and store them into a hash table
   -----------------------------------------------------------------
*/
   ALLOCATE(frontmtx->p_mtxDJJ, struct _SubMtx *, nfront) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      frontmtx->p_mtxDJJ[J] = NULL ;
   }
   if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n nmtx = %d", nmtx) ;
   fflush(stdout) ;
#endif
   for ( J = 0 ; J < nmtx ; J++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromBinaryFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in diag SubMtx object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
      frontmtx->p_mtxDJJ[mtx->rowid] = mtx ;
#if MYDEBUG > 0
      fprintf(stdout, "\n\n --- diagonal SubMtx object") ;
      SubMtx_writeForHumanEye(mtx, stdout) ;
      fflush(stdout) ;
#endif
   }
/*
   --------------------------------------------------------------
   read in the upper submatrices and store them into a hash table
   --------------------------------------------------------------
*/
   frontmtx->upperhash = I2Ohash_new() ;
   if ( (rc = fread((void *) &nmtx, sizeof(int), 1, fp)) != 1 ) {
      fprintf(stderr, "\n error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n %d items of %d read\n", frontmtx, fp, rc, 1) ;
      return(0) ;
   }
#if MYDEBUG > 0
   fprintf(stdout, "\n\n nmtx = %d", nmtx) ;
   fflush(stdout) ;
#endif
   I2Ohash_init(frontmtx->upperhash, nfront-1, nmtx, 0) ;
   for ( JK = 0 ; JK < nmtx ; JK++ ) {
      mtx = SubMtx_new() ;
      rc = SubMtx_readFromBinaryFile(mtx, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr,
              "\n fatal error in FrontMtx_readFromBinaryFile(%p,%p)"
              "\n error %d reading in upper SubMtx object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
      I2Ohash_insert(frontmtx->upperhash, 
                     mtx->rowid, mtx->colid, (void *) mtx) ;
#if MYDEBUG > 0
      fprintf(stdout, "\n\n --- upper SubMtx object") ;
      SubMtx_writeForHumanEye(mtx, stdout) ;
      fflush(stdout) ;
#endif
   }
}
#if MYDEBUG > 0
   fprintf(stdout, "\n\n LEAVING READ_FROM_BINARY_FILE") ;
   fflush(stdout) ;
#endif
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to write an FrontMtx object to a file

   input --

      fn -- filename
        *.frontmtxb -- binary
        *.frontmtxf -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 98may04, cca
   -------------------------------------------------
*/
int
FrontMtx_writeToFile ( 
   FrontMtx   *frontmtx, 
   char        *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_writeToFile(%p,%s)"
    "\n bad input\n", frontmtx, fn) ; 
}
/*
   ------------------
   write out the file
   ------------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "wb")) == NULL ) {
         fprintf(stderr, "\n error in FrontMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", frontmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = FrontMtx_writeToBinaryFile(frontmtx, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in FrontMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", frontmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = FrontMtx_writeToFormattedFile(frontmtx, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in FrontMtx_writeToFile(%p,%s)"
                 "\n unable to open file %s", frontmtx, fn, fn) ;
         rc = 0 ;
      } else {
         rc = FrontMtx_writeForHumanEye(frontmtx, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in FrontMtx_writeToFile(%p,%s)"
              "\n unable to open file %s", frontmtx, fn, fn) ;
      rc = 0 ;
   } else {
      rc = FrontMtx_writeForHumanEye(frontmtx, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write an FrontMtx object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 98may04, cca
   -----------------------------------------------------------
*/
int
FrontMtx_writeToFormattedFile ( 
   FrontMtx   *frontmtx, 
   FILE       *fp 
) {
SubMtx   *mtx ;
int    ii, J, K, nadj, nfront, nmtx, rc ;
int    *adj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
/*
   ---------------------------
   write the ten scalar fields
   ---------------------------
*/
rc = fprintf(fp, "\n %d %d %d %d %d %d %d %d %d %d", 
             frontmtx->nfront, frontmtx->neqns, frontmtx->type,
             frontmtx->symmetryflag, frontmtx->pivotingflag, 
             frontmtx->sparsityflag, frontmtx->dataMode,
             frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", 
           frontmtx, fp, rc) ;
   return(0) ;
}
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
rc = ETree_writeToFormattedFile(frontmtx->frontETree, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
           "\n error %d writing frontETree object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
/*
   -------------------------------------------
   write out the symbolic factorization object
   -------------------------------------------
*/
rc = IVL_writeToFormattedFile(frontmtx->symbfacIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
           "\n error %d writing symbfacIVL object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   --------------------------------
   write out the front sizes object
   --------------------------------
*/
   rc = IV_writeToFormattedFile(frontmtx->frontsizesIV, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
              "\n error %d writing frontsizesIV object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   -------------------------------------------
   write out the rowids and colids IVL objects
   -------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      rc = IVL_writeToFormattedFile(frontmtx->rowadjIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
              "\n error %d writing rowadjIVL object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   rc = IVL_writeToFormattedFile(frontmtx->coladjIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
           "\n error %d writing coladjIVL object\n",
           frontmtx, fp, rc) ;
      return(0) ;
   }
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      ------------------------------------------------
      count the number of L_{J,J} matrices
      write out the count, then write out the matrices
      ------------------------------------------------
*/
      for ( J = nmtx = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, J, J)) != NULL ) {
            nmtx++ ;
         }
      }
      fprintf(fp, "\n %d", nmtx) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, J, J)) != NULL ) {
            SubMtx_writeToFormattedFile(mtx, fp) ;
         }
      }
/*
      ------------------------------------------------
      count the number of L_{*,J} matrices
      write out the count, then write out the matrices
      ------------------------------------------------
*/
      for ( J = nmtx = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, nfront, J)) != NULL ){
            nmtx++ ;
         }
      }
      fprintf(fp, "\n %d", nmtx) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, nfront, J)) != NULL ){
            SubMtx_writeToFormattedFile(mtx, fp) ;
         }
      }
   }
/*
   ------------------------------------------------
   count the number of diagonal matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         nmtx++ ;
      }
   }
   fprintf(fp, "\n %d", nmtx) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         SubMtx_writeToFormattedFile(mtx, fp) ;
      }
   }
/*
   ------------------------------------------------
   count the number of U_{J,J} matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
         nmtx++ ;
      }
   }
   fprintf(fp, "\n %d", nmtx) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
         SubMtx_writeToFormattedFile(mtx, fp) ;
      }
   }
/*
   ------------------------------------------------
   count the number of U_{J,*} matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, nfront)) != NULL ) {
         nmtx++ ;
      }
   }
   fprintf(fp, "\n %d", nmtx) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, nfront)) != NULL ) {
         SubMtx_writeToFormattedFile(mtx, fp) ;
      }
   }
} else {
/*
   ---------------------------------------------------------
   write out the lower and upper block adjacency IVL objects
   ---------------------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      rc = IVL_writeToFormattedFile(frontmtx->lowerblockIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
              "\n error %d writing lowerblockIVL object\n",
              frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   rc = IVL_writeToFormattedFile(frontmtx->upperblockIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToFormattedFile(%p,%p)"
              "\n error %d writing upperblockIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -------------------------------
      write out the lower submatrices
      -------------------------------
*/
      for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
         IVL_listAndSize(frontmtx->lowerblockIVL, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            K = adj[ii] ;
            if ( (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               nmtx++ ;
            }
         }
      }
      fprintf(fp, "\n %d", nmtx) ;
      for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
         IVL_listAndSize(frontmtx->lowerblockIVL, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            K = adj[ii] ;
            if ( (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               SubMtx_writeToFormattedFile(mtx, fp) ;
            }
         }
      }
   }
/*
   ----------------------------------
   write out the diagonal submatrices
   ----------------------------------
*/
   for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         nmtx++ ;
      }
   }
   fprintf(fp, "\n %d", nmtx) ;
   for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         SubMtx_writeToFormattedFile(mtx, fp) ;
      }
   }
/*
   -------------------------------
   write out the upper submatrices
   -------------------------------
*/
   for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
      IVL_listAndSize(frontmtx->upperblockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = adj[ii] ;
         if ( (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            nmtx++ ;
         }
      }
   }
   fprintf(fp, "\n %d", nmtx) ;
   for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
      IVL_listAndSize(frontmtx->upperblockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = adj[ii] ;
         if ( (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            SubMtx_writeToFormattedFile(mtx, fp) ;
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   purpose -- to write an FrontMtx object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 98may04, cca
   -------------------------------------------------------
*/
int
FrontMtx_writeToBinaryFile ( 
   FrontMtx   *frontmtx, 
   FILE        *fp 
) {
SubMtx   *mtx ;
int    ii, J, K, nadj, nfront, nmtx, rc ;
int    *adj ;
int    itemp[10] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   exit(-1) ;
}
/*
   ---------------------------
   write the ten scalar fields
   ---------------------------
*/
itemp[0] = nfront = frontmtx->nfront ;
itemp[1] = frontmtx->neqns ;
itemp[2] = frontmtx->type ;
itemp[3] = frontmtx->symmetryflag ;
itemp[4] = frontmtx->pivotingflag ;
itemp[5] = frontmtx->sparsityflag ;
itemp[6] = frontmtx->dataMode  ;
itemp[7] = frontmtx->nentD  ;
itemp[8] = frontmtx->nentL  ;
itemp[9] = frontmtx->nentU  ;
rc = fwrite((void *) itemp, sizeof(int), 10, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", 
           frontmtx, fp, rc) ;
   return(0) ;
}
/*
   --------------------------
   write out the ETree object
   --------------------------
*/
rc = ETree_writeToBinaryFile(frontmtx->frontETree, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
           "\n error %d writing frontETree object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
/*
   -------------------------------------------
   write out the symbolic factorization object
   -------------------------------------------
*/
rc = IVL_writeToBinaryFile(frontmtx->symbfacIVL, fp) ;
if ( rc != 1 ) {
   fprintf(stderr, 
           "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
           "\n error %d writing symbfacIVL object\n",
           frontmtx, fp, rc) ;
   return(0) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   --------------------------------
   write out the front sizes object
   --------------------------------
*/
   rc = IV_writeToBinaryFile(frontmtx->frontsizesIV, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
              "\n error %d writing frontsizesIV object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
/*
   -------------------------------------------
   write out the rowids and colids IVL objects
   -------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      rc = IVL_writeToBinaryFile(frontmtx->rowadjIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
                 "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
                 "\n error %d writing rowadjIVL object\n",
                 frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   rc = IVL_writeToBinaryFile(frontmtx->coladjIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
              "\n error %d writing coladjIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      ------------------------------------------------
      count the number of L_{J,J} matrices
      write out the count, then write out the matrices
      ------------------------------------------------
*/
      for ( J = nmtx = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, J, J)) != NULL ) {
            nmtx++ ;
         }
      }
      rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, J, J)) != NULL ) {
            SubMtx_writeToBinaryFile(mtx, fp) ;
         }
      }
/*
      ------------------------------------------------
      count the number of L_{*,J} matrices
      write out the count, then write out the matrices
      ------------------------------------------------
*/
      for ( J = nmtx = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, nfront, J)) != NULL ){
            nmtx++ ;
         }
      }
      rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
      for ( J = 0 ; J < nfront ; J++ ) {
         if ( (mtx = FrontMtx_lowerMtx(frontmtx, nfront, J)) != NULL ){
            SubMtx_writeToBinaryFile(mtx, fp) ;
         }
      }
   }
/*
   ------------------------------------------------
   count the number of diagonal matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         nmtx++ ;
      }
   }
   rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         SubMtx_writeToBinaryFile(mtx, fp) ;
      }
   }
/*
   ------------------------------------------------
   count the number of U_{J,J} matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
         nmtx++ ;
      }
   }
   rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, J)) != NULL ) {
         SubMtx_writeToBinaryFile(mtx, fp) ;
      }
   }
/*
   ------------------------------------------------
   count the number of U_{J,*} matrices
   write out the count, then write out the matrices
   ------------------------------------------------
*/
   for ( J = nmtx = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, nfront)) != NULL ) {
         nmtx++ ;
      }
   }
   rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
   for ( J = 0 ; J < nfront ; J++ ) {
      if ( (mtx = FrontMtx_upperMtx(frontmtx, J, nfront)) != NULL ) {
         SubMtx_writeToBinaryFile(mtx, fp) ;
      }
   }
} else {
/*
   ---------------------------------------------------------
   write out the lower and upper block adjacency IVL objects
   ---------------------------------------------------------
*/
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      rc = IVL_writeToBinaryFile(frontmtx->lowerblockIVL, fp) ;
      if ( rc != 1 ) {
         fprintf(stderr, 
                 "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
                 "\n error %d writing lowerblockIVL object\n",
                 frontmtx, fp, rc) ;
         return(0) ;
      }
   }
   rc = IVL_writeToBinaryFile(frontmtx->upperblockIVL, fp) ;
   if ( rc != 1 ) {
      fprintf(stderr, 
              "\n fatal error in FrontMtx_writeToBinaryFile(%p,%p)"
              "\n error %d writing upperblockIVL object\n",
              frontmtx, fp, rc) ;
      return(0) ;
   }
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
/*
      -------------------------------
      write out the lower submatrices
      -------------------------------
*/
      for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
         IVL_listAndSize(frontmtx->lowerblockIVL, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            K = adj[ii] ;
            if ( (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               nmtx++ ;
            }
         }
      }
      rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
      for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
         IVL_listAndSize(frontmtx->lowerblockIVL, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            K = adj[ii] ;
            if ( (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               SubMtx_writeToBinaryFile(mtx, fp) ;
            }
         }
      }
   }
/*
   ----------------------------------
   write out the diagonal submatrices
   ----------------------------------
*/
   for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         nmtx++ ;
      }
   }
   rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
   for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
      if ( (mtx = FrontMtx_diagMtx(frontmtx, J)) != NULL ) {
         SubMtx_writeToBinaryFile(mtx, fp) ;
      }
   }
/*
   -------------------------------
   write out the upper submatrices
   -------------------------------
*/
   for ( J = nmtx = 0 ; J < frontmtx->nfront ; J++ ) {
      IVL_listAndSize(frontmtx->upperblockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = adj[ii] ;
         if ( (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            nmtx++ ;
         }
      }
   }
   rc = fwrite((void *) &nmtx, sizeof(int), 1, fp) ;
   for ( J = 0 ; J < frontmtx->nfront ; J++ ) {
      IVL_listAndSize(frontmtx->upperblockIVL, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         K = adj[ii] ;
         if ( (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            SubMtx_writeToBinaryFile(mtx, fp) ;
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   purpose -- to write out the statistics for the FrontMtx object

   return value -- 1 if success, 0 otherwise

   created -- 98may04, cca
   ---------------------------------------------------------------
*/
int
FrontMtx_writeStats ( 
   FrontMtx   *frontmtx, 
   FILE        *fp 
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in FrontMtx_writeStats(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   exit(-1) ;
}
fprintf(fp, "\n\n FrontMtx object at address %p"
        "\n nfront = %d", frontmtx, frontmtx->nfront) ;
switch ( frontmtx->symmetryflag ) {
case SPOOLES_SYMMETRIC : 
   fprintf(fp, "\n symmetric entries") ;
   break ;
case SPOOLES_HERMITIAN : 
   fprintf(fp, "\n Hermitian") ;
   break ;
case SPOOLES_NONSYMMETRIC : 
   fprintf(fp, "\n nonsymmetric structure, nonsymmetric entries") ;
   break ;
default :
   break ;
}
switch ( frontmtx->pivotingflag ) {
case SPOOLES_NO_PIVOTING : 
   fprintf(fp, "\n pivoting disabled") ;
   break ;
case SPOOLES_PIVOTING : 
   fprintf(fp, "\n pivoting enabled") ;
   break ;
default :
   break ;
}
switch ( frontmtx->sparsityflag ) {
case FRONTMTX_DENSE_FRONTS : 
   fprintf(fp, "\n dense fronts") ;
   break ;
case FRONTMTX_SPARSE_FRONTS : 
   fprintf(fp, "\n sparse fronts ") ;
   break ;
default :
   break ;
}
switch ( frontmtx->dataMode ) {
case FRONTMTX_1D_MODE : 
   fprintf(fp, "\n one-dimensional data decomposition") ;
   break ;
case FRONTMTX_2D_MODE : 
   fprintf(fp, "\n two-dimensional data decomposition") ;
   break ;
default :
   break ;
}
fprintf(fp, "\n %d entries in D, %d entries in L, %d entries in U",
        frontmtx->nentD, frontmtx->nentL, frontmtx->nentU) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98may04, cca
   ----------------------------------------
*/
int
FrontMtx_writeForHumanEye (
   FrontMtx   *frontmtx,
   FILE       *fp
) {
SubMtx   *mtx ;
int    ii, J, K, nadj, nfront ;
int    *adj ;
/*
   ---------------
   check the input
   ---------------
*/
if ( frontmtx == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_writeForHumanEye(%p,%p)"
           "\n bad input\n", frontmtx, fp) ;
   exit(-1) ;
}
nfront = frontmtx->nfront ;
FrontMtx_writeStats(frontmtx, fp) ;
if ( frontmtx->frontETree != NULL ) {
   fprintf(fp, "\n\n front tree FrontMtx object") ;
   ETree_writeForHumanEye(frontmtx->frontETree, fp) ;
   fflush(fp) ;
}
if ( frontmtx->symbfacIVL != NULL ) {
   fprintf(fp, "\n\n symbfacIVL IVL object") ;
   IVL_writeForHumanEye(frontmtx->symbfacIVL, fp) ;
   fflush(fp) ;
}
if ( frontmtx->frontsizesIV != NULL ) {
   fprintf(fp, "\n\n frontsizesIV IV object") ;
   IV_writeForHumanEye(frontmtx->frontsizesIV, fp) ;
   fflush(fp) ;
}
if ( frontmtx->rowadjIVL != NULL ) {
   fprintf(fp, "\n\n rowids IVL object") ;
   IVL_writeForHumanEye(frontmtx->rowadjIVL, fp) ;
   fflush(fp) ;
}
if ( frontmtx->coladjIVL != NULL ) {
   fprintf(fp, "\n\n colids IVL object") ;
   IVL_writeForHumanEye(frontmtx->coladjIVL, fp) ;
   fflush(fp) ;
}
if ( frontmtx->lowerblockIVL != NULL ) {
   fprintf(fp, "\n\n lower block adjacency IVL object") ;
   IVL_writeForHumanEye(frontmtx->lowerblockIVL, fp) ;
   fflush(fp) ;
}
if ( frontmtx->upperblockIVL != NULL ) {
   fprintf(fp, "\n\n upper block adjacency IVL object") ;
   IVL_writeForHumanEye(frontmtx->upperblockIVL, fp) ;
   fflush(fp) ;
}
if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
   fprintf(fp, "\n\n lower submatrices") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      mtx = FrontMtx_lowerMtx(frontmtx, J, J) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n --- lower submatrix -- diagonal") ;
         SubMtx_writeForHumanEye(mtx, fp) ;
         fflush(fp) ;
      }
      if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
         mtx = FrontMtx_lowerMtx(frontmtx, nfront, J) ;
         if ( mtx != NULL ) {
            fprintf(fp, "\n\n --- lower submatrix") ;
            SubMtx_writeForHumanEye(mtx, fp) ;
            fflush(fp) ;
         }
      } else {
         FrontMtx_lowerAdjFronts(frontmtx, J, &nadj, &adj) ;
         for ( ii = 0 ; ii < nadj ; ii++ ) {
            if (  (K = adj[ii]) > J 
               && (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               fprintf(fp, "\n\n --- lower submatrix") ;
               SubMtx_writeForHumanEye(mtx, fp) ;
               fflush(fp) ;
            }
         }
      }
   }
}
fprintf(fp, "\n\n diagonal submatrices") ;
for ( J = 0 ; J < nfront ; J++ ) {
   mtx = FrontMtx_diagMtx(frontmtx, J) ;
   if ( mtx != NULL ) {
      fprintf(fp, "\n\n --- diagonal submatrix") ;
      SubMtx_writeForHumanEye(mtx, fp) ;
   }
   fflush(fp) ;
}
fprintf(fp, "\n\n upper submatrices") ;
for ( J = 0 ; J < nfront ; J++ ) {
   mtx = FrontMtx_upperMtx(frontmtx, J, J) ;
   if ( mtx != NULL ) {
      fprintf(fp, "\n\n --- upper submatrix --- diagonal") ;
      SubMtx_writeForHumanEye(mtx, fp) ;
      fflush(fp) ;
   }
   if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
      mtx = FrontMtx_upperMtx(frontmtx, J, nfront) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n --- upper submatrix") ;
         SubMtx_writeForHumanEye(mtx, fp) ;
         fflush(fp) ;
      }
   } else {
      FrontMtx_upperAdjFronts(frontmtx, J, &nadj, &adj) ;
      for ( ii = 0 ; ii < nadj ; ii++ ) {
         if (  (K = adj[ii]) > J 
            && (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            fprintf(fp, "\n\n --- upper submatrix") ;
            SubMtx_writeForHumanEye(mtx, fp) ;
            fflush(fp) ;
         }
      }
   }
}
fprintf(fp, "\n\n ### leaving FrontMtx_writeForHumanEye") ;
fflush(fp) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------------
   purpose -- to write the factor matrices out for a matlab file

      Lname -- name for lower triangular matrix
      Dname -- name for diagonal matrix
      Uname -- name for upper triangular matrix

   presently works only with 1-d data decomposition

   created -- 98sep23, cca
   -------------------------------------------------------------
*/
int 
FrontMtx_writeForMatlab (
   FrontMtx   *frontmtx,
   char       *Lname,
   char       *Dname,
   char       *Uname,
   FILE       *fp
) {
int      J, nfront ;
SubMtx   *mtx ;
/*
   ---------------
   check the input
   ---------------
*/
if (  frontmtx == NULL || Lname == NULL || Dname == NULL 
   || Uname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in FrontMtx_writeForMatlab()"
           "\n bad input\n") ;
   exit(-1) ;
}
if ( FRONTMTX_IS_1D_MODE(frontmtx) ) {
   nfront = FrontMtx_nfront(frontmtx) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      fprintf(fp, "\n\n %% lower submatrices") ;
      for ( J = 0 ; J < nfront ; J++ ) {
         mtx = FrontMtx_lowerMtx(frontmtx, J, J) ;
         if ( mtx != NULL ) {
            fprintf(fp, "\n\n %% --- lower submatrix -- diagonal") ;
            SubMtx_writeForMatlab(mtx, Lname, fp) ;
            fflush(fp) ;
         }
         mtx = FrontMtx_lowerMtx(frontmtx, nfront, J) ;
         if ( mtx != NULL ) {
            fprintf(fp, "\n\n %% --- lower submatrix") ;
            SubMtx_writeForMatlab(mtx, Lname, fp) ;
            fflush(fp) ;
         }
      }
   }
   fprintf(fp, "\n\n %% diagonal submatrices") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      mtx = FrontMtx_diagMtx(frontmtx, J) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n %% --- diagonal submatrix") ;
         SubMtx_writeForMatlab(mtx, Dname, fp) ;
      }
      fflush(fp) ;
   }
   fprintf(fp, "\n\n %% upper submatrices") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      mtx = FrontMtx_upperMtx(frontmtx, J, J) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n %% --- upper submatrix --- diagonal") ;
         SubMtx_writeForMatlab(mtx, Uname, fp) ;
         fflush(fp) ;
      }
      mtx = FrontMtx_upperMtx(frontmtx, J, nfront) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n %% --- upper submatrix") ;
         SubMtx_writeForMatlab(mtx, Uname, fp) ;
         fflush(fp) ;
      }
   }
} else if ( FRONTMTX_IS_2D_MODE(frontmtx) ) {
   nfront = FrontMtx_nfront(frontmtx) ;
   if ( FRONTMTX_IS_NONSYMMETRIC(frontmtx) ) {
      int   ii, jj, kk, K, nadj, ncolJ, ncolKJ, nrowK, nrowKJ ;
      int   *adj, *colindJ, *colKJ, *rowindK, *rowKJ ;

      fprintf(fp, "\n\n %% lower submatrices") ;
      for ( J = 0 ; J < nfront ; J++ ) {
         mtx = FrontMtx_lowerMtx(frontmtx, J, J) ;
         if ( mtx != NULL ) {
            fprintf(fp, "\n\n %% --- lower submatrix -- diagonal") ;
            SubMtx_writeForMatlab(mtx, Lname, fp) ;
            fflush(fp) ;
         }
         FrontMtx_lowerAdjFronts(frontmtx, J, &nadj, &adj) ;
         for ( kk = 0 ; kk < nadj ; kk++ ) {
            if (  (K = adj[kk]) > J 
               && (mtx = FrontMtx_lowerMtx(frontmtx, K, J)) != NULL ) {
               fprintf(fp, "\n\n %% --- lower submatrix") ;
               SubMtx_columnIndices(mtx, &ncolKJ, &colKJ) ;
               FrontMtx_columnIndices(frontmtx, J, &ncolJ, &colindJ) ;
               for ( ii = 0 ; ii < ncolKJ ; ii++ ) {
                  colKJ[ii] = colindJ[colKJ[ii]] ;
               }
               SubMtx_rowIndices(mtx, &nrowKJ, &rowKJ) ;
               FrontMtx_rowIndices(frontmtx, K, &nrowK, &rowindK) ;
               for ( ii = 0 ; ii < nrowKJ ; ii++ ) {
                  rowKJ[ii] = rowindK[rowKJ[ii]] ;
               }
               SubMtx_writeForMatlab(mtx, Lname, fp) ;
               for ( ii = jj = 0 ; ii < ncolKJ ; ii++ ) {
                  while ( colKJ[ii] != colindJ[jj] ) {
                     jj++ ;
                  }
                  colKJ[ii] = jj++ ;
               }
               for ( ii = jj = 0 ; ii < nrowKJ ; ii++ ) {
                  while ( rowKJ[ii] != rowindK[jj] ) {
                     jj++ ;
                  }
                  rowKJ[ii] = jj++ ;
               }
               fflush(fp) ;
            }
         }
      }
   }
   fprintf(fp, "\n\n %% diagonal submatrices") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      mtx = FrontMtx_diagMtx(frontmtx, J) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n %% --- diagonal submatrix") ;
         SubMtx_writeForMatlab(mtx, Dname, fp) ;
      }
      fflush(fp) ;
   }
   fprintf(fp, "\n\n %% upper submatrices") ;
   for ( J = 0 ; J < nfront ; J++ ) {
      int   ii, jj, kk, K, nadj, ncolK, ncolJK, nrowJ, nrowJK ;
      int   *adj, *colindK, *colJK, *rowindJ, *rowJK ;

      mtx = FrontMtx_upperMtx(frontmtx, J, J) ;
      if ( mtx != NULL ) {
         fprintf(fp, "\n\n %% --- upper submatrix --- diagonal") ;
         SubMtx_writeForMatlab(mtx, Uname, fp) ;
         fflush(fp) ;
      }
      FrontMtx_upperAdjFronts(frontmtx, J, &nadj, &adj) ;
      for ( kk = 0 ; kk < nadj ; kk++ ) {
         if (  (K = adj[kk]) > J 
            && (mtx = FrontMtx_upperMtx(frontmtx, J, K)) != NULL ) {
            fprintf(fp, "\n\n %% --- upper submatrix") ;
            SubMtx_columnIndices(mtx, &ncolJK, &colJK) ;
            FrontMtx_columnIndices(frontmtx, K, &ncolK, &colindK) ;
            for ( ii = 0 ; ii < ncolJK ; ii++ ) {
               colJK[ii] = colindK[colJK[ii]] ;
            }
            SubMtx_rowIndices(mtx, &nrowJK, &rowJK) ;
            FrontMtx_rowIndices(frontmtx, J, &nrowJ, &rowindJ) ;
            for ( ii = 0 ; ii < nrowJK ; ii++ ) {
               rowJK[ii] = rowindJ[rowJK[ii]] ;
            }
            SubMtx_writeForMatlab(mtx, Uname, fp) ;
            for ( ii = jj = 0 ; ii < ncolJK ; ii++ ) {
               while ( colJK[ii] != colindK[jj] ) {
                  jj++ ;
               }
               colJK[ii] = jj++ ;
            }
            for ( ii = jj = 0 ; ii < nrowJK ; ii++ ) {
               while ( rowJK[ii] != rowindJ[jj] ) {
                  jj++ ;
               }
               rowJK[ii] = jj++ ;
            }
            fflush(fp) ;
         }
      }
   }
}
return(1) ; }

/*--------------------------------------------------------------------*/
