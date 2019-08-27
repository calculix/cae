/*  instance.c  */

#include "../Chv.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   return the id of the chevron

   created -- 98apr30, cca
   ----------------------------
*/
int
Chv_id (
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_id(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return(chv->id) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   return the type of the chevron
   return value = SPOOLES_REAL    --> chevron is real
   return value = SPOOLES_COMPLEX --> chevron is complex

   created -- 98apr30, cca
   ------------------------------------------------------------
*/
int
Chv_type (
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_type(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return(chv->type) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------------
   return the symmetry flag of the chevron
   return value = SPOOLES_SYMMETRIC --> chevron is symmetric
   return value = SPOOLES_HERMITIAN --> chevron is hermitian
   return value = SPOOLES_NONSYMMETRIC --> chevron is nonsymmetric

   created -- 98apr30, cca
   ------------------------------------------------------------
*/
int
Chv_symmetryFlag (
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_symmetryFlag(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return(chv->symflag) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   fill *pnD with nD, *pnL with nL, and *pnU with nU.

   created -- 98apr30, cca
   --------------------------------------------------
*/
void
Chv_dimensions (
   Chv   *chv,
   int   *pnD,
   int   *pnL,
   int   *pnU
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || pnD == NULL || pnL == NULL || pnU == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_dimensions(%p,%p,%p,%p)"
           "\n bad input\n", chv, pnD, pnL, pnU) ;
   exit(-1) ;
}
*pnD = chv->nD ;
*pnL = chv->nL ;
*pnU = chv->nU ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   fill *pnrow with nD + nL, *prowind with rowind

   created -- 98apr30, cca
   ----------------------------------------------
*/
void
Chv_rowIndices (
   Chv   *chv,
   int   *pnrow,
   int   **prowind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || pnrow == NULL || prowind == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_rowIndices(%p,%p,%p)"
           "\n bad input\n", chv, pnrow, prowind) ;
   exit(-1) ;
}
if ( CHV_IS_NONSYMMETRIC(chv) ) {
   *pnrow   = chv->nD + chv->nL ;
   *prowind = chv->rowind ;
} else if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   *pnrow   = chv->nD + chv->nU ;
   *prowind = chv->colind ;
} else {
   fprintf(stderr, "\n fatal error in Chv_rowIndices(%p,%p,%p)"
           "\n bad symflag = %d\n", chv, pnrow, prowind, chv->symflag) ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   fill *pncol with nD + nU, *pcolind with colind

   created -- 98apr30, cca
   ----------------------------------------------
*/
void
Chv_columnIndices (
   Chv   *chv,
   int   *pncol,
   int   **pcolind
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || pncol == NULL || pcolind == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_columnIndices(%p,%p,%p)"
           "\n bad input\n", chv, pncol, pcolind) ;
   exit(-1) ;
}
*pncol   = chv->nD + chv->nU ;
*pcolind = chv->colind ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------
   return the number of entries

   created -- 98apr30, cca
   ----------------------------
*/
int  
Chv_nent (
   Chv   *chv
) {
int   nD, nent, nL, nU ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_nent(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   nent = (nD*(nD+1))/2 + nD*nU ;
} else if ( CHV_IS_NONSYMMETRIC(chv) ) {
   nent = nD*(nD + nL + nU) ;
} else {
   fprintf(stderr, "\n fatal error in Chv_nent(%p)"
           "\n bad symmetry flag %d\n", chv, chv->symflag) ;
   exit(-1) ;
}
return(nent) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill *pentries with a pointer to the entries

   created -- 98apr30, cca
   --------------------------------------------
*/
double *
Chv_entries(
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_entries(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return(chv->entries) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   return the location of the diagonal entry
   for the ichv'th chevron

   created -- 98apr30, cca
   -----------------------------------------
*/
double *
Chv_diagLocation(
   Chv   *chv,
   int    ichv
) {
double   *diag ;
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
if ( ichv < 0 || ichv > chv->nD ) {
   fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
           "\n ichv = %d, nD = %d\n", chv, ichv, chv->nD) ;
   exit(-1) ;
}
if ( chv->entries == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
           "\n chv->entries is NULL\n", chv) ;
   exit(-1) ;
}
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
      diag = chv->entries + ichv*(chv->nD + chv->nU) 
                                  - (ichv*(ichv-1))/2 ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      diag = chv->entries + (2*ichv+1)*chv->nD + (ichv+1)*chv->nL
           + ichv*chv->nU - ichv*ichv - ichv - 1 ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
              "\n type is SPOOLES_REAL, symflag = %d" 
              "\n not SPOOLES_SYMMETRIC or SPOOLES_NONSYMMETRIC\n", 
              chv, chv->symflag) ;
      exit(-1) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
      diag = chv->entries + 2*(ichv*(chv->nD + chv->nU) 
                                  - (ichv*(ichv-1))/2) ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      diag = chv->entries + 2*((2*ichv+1)*chv->nD + (ichv+1)*chv->nL
           + ichv*chv->nU - ichv*ichv - ichv - 1) ;
   } else {
      fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
              "\n bad symflag = %d, type is SPOOLES_COMPLEX,"
              "\n must be SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
              "\n or SPOOLES_NONSYMMETRIC\n",
              chv, chv->symflag) ;
      exit(-1) ;
   }
} else {
   fprintf(stderr, "\n fatal error in Chv_diagLocation(%p)"
           "\n bad type = %d, not SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, chv->symflag) ;
   exit(-1) ;
}
return(diag) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return a pointer to the start of the workspace

   created -- 98apr30, cca
   ----------------------------------------------
*/
void *
Chv_workspace(
   Chv   *chv
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  chv == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_workspace(%p)"
           "\n bad input\n", chv) ;
   exit(-1) ;
}
return((void *) DV_entries(&chv->wrkDV)) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   fill *pValue with entry (irow, jcol)
 
   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_realEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   *pValue
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 
    || pValue == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_realEntry(%p,%d,%d,%p)"
           "\n bad input\n", chv, irow, jcol, pValue) ;
   exit(-1) ;
}
if ( ! CHV_IS_REAL(chv) ) {
   fprintf(stderr, "\n fatal error in Chv_realEntry(%p,%d,%d,%p)"
           "\n bad type %d, not SPOOLES_REAL\n", 
           chv, irow, jcol, pValue, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_realEntry(%p,%d,%d,%p)"
           "\n bad symflag %d"
           "\n must be SPOOLES_SYMMETRIC of SPOOLES_NONSYMMETRIC\n",
           chv, irow, jcol, pValue, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, "\n fatal error in Chv_realEntry(%p,%d,%d,%p)"
           "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
           chv, irow, jcol, pValue, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow >= nD && jcol >= nD ) {
   *pValue = 0.0 ;
} else {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( CHV_IS_SYMMETRIC(chv) && off < 0 ) {
      off = -off ;
   }
   base    = Chv_diagLocation(chv, ichv) ;
   *pValue = base[off] ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill (*pReal,*pImag) with entry (irow, jcol)
 
   created -- 98apr30, cca
   --------------------------------------------
*/
void
Chv_complexEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   *pReal,
   double   *pImag
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 
    || pReal == NULL || pImag == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad input\n", chv, irow, jcol, pReal, pImag) ;
   exit(-1) ;
}
if ( ! CHV_IS_COMPLEX(chv) ) {
   fprintf(stderr, "\n fatal error in Chv_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad type %d, not SPOOLES_COMPLEX\n", 
           chv, irow, jcol, pReal, pImag, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) 
        || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_complexEntry(%p,%d,%d,%p,%p)"
           "\n bad symflag %d, not SPOOLES_SYMMETRIC, "
           "\n SPOOLES_HERMITIAN or SPOOLES_NONSYMMETRIC \n",
           chv, irow, jcol, pReal, pImag, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, "\n fatal error in Chv_complexEntry(%p,%d,%d,%p,%p)"
           "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
           chv, irow, jcol, pReal, pImag, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow >= nD && jcol >= nD ) {
   *pReal = *pImag = 0.0 ;
} else {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( off < 0 && (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv)) ) {
      off = -off ;
   }
   base   = Chv_diagLocation(chv, ichv) ;
   *pReal = base[2*off]   ;
   if ( irow > jcol && CHV_IS_HERMITIAN(chv) ) {
      *pImag = - base[2*off+1] ;
   } else {
      *pImag = base[2*off+1] ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------
   fill *ppValue with the location of entry (irow, jcol)
 
   created -- 98apr30, cca
   -----------------------------------------------------
*/
void
Chv_locationOfRealEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   **ppValue
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 
    || ppValue == NULL ) {
   fprintf(stderr, 
           "\n fatal error in Chv_locationOfRealEntry(%p,%d,%d,%p)"
           "\n bad input\n", chv, irow, jcol, ppValue) ;
   exit(-1) ;
}
if ( ! CHV_IS_REAL(chv) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_locationOfRealEntry(%p,%d,%d,%p)"
           "\n bad type %d, not SPOOLES_REAL\n", 
           chv, irow, jcol, ppValue, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_locationOfRealEntry(%p,%d,%d,%p)"
           "\n bad symflag %d, not SPOOLES_SYMMETRIC of SPOOLES_NONSYMMETRIC\n",
           chv, irow, jcol, ppValue, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, 
           "\n fatal error in Chv_locationOfRealEntry(%p,%d,%d,%p)"
           "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
           chv, irow, jcol, ppValue, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow >= nD && jcol >= nD ) {
   *ppValue = NULL ;
} else {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( CHV_IS_SYMMETRIC(chv) && off < 0 ) {
      off = -off ;
   }
   base     = Chv_diagLocation(chv, ichv) ;
   *ppValue = base + off ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------------
   fill (*ppReal,*ppImag) with location of entry (irow, jcol)
 
   created -- 98apr30, cca
   ----------------------------------------------------------
*/
void
Chv_locationOfComplexEntry (
   Chv     *chv,
   int      irow,
   int      jcol,
   double   **ppReal,
   double   **ppImag
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 
    || ppReal == NULL || ppImag == NULL ) {
   fprintf(stderr, 
          "\n fatal error in Chv_locationOfComplexEntry(%p,%d,%d,%p,%p)"
          "\n bad input\n", chv, irow, jcol, ppReal, ppImag) ;
   exit(-1) ;
}
if ( ! CHV_IS_COMPLEX(chv) ) {
   fprintf(stderr, 
          "\n fatal error in Chv_locationOfComplexEntry(%p,%d,%d,%p,%p)"
          "\n bad type %d, not SPOOLES_COMPLEX\n", 
          chv, irow, jcol, ppReal, ppImag, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) 
        || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, 
          "\n fatal error in Chv_locationOfComplexEntry(%p,%d,%d,%p,%p)"
          "\n bad symflag %d"
          "\n not SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
          "\n or SPOOLES_NONSYMMETRIC \n",
          chv, irow, jcol, ppReal, ppImag, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, 
          "\n fatal error in Chv_locationOfComplexEntry(%p,%d,%d,%p,%p)"
          "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
          chv, irow, jcol, ppReal, ppImag, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow >= nD && jcol >= nD ) {
   *ppReal = *ppImag = NULL ;
} else {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( off < 0 && (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv)) ) {
      off = -off ;
   }
   base    = Chv_diagLocation(chv, ichv) ;
   *ppReal = base + 2*off ;
   *ppImag = base + 2*off + 1 ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------
   set entry (irow, jcol) to value
 
   created -- 98apr30, cca
   ------------------------------------
*/
void
Chv_setRealEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   value
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 ) {
   fprintf(stderr, "\n fatal error in Chv_setRealEntry(%p,%d,%d,%e)"
           "\n bad input\n", chv, irow, jcol, value) ;
   exit(-1) ;
}
if ( ! CHV_IS_REAL(chv) ) {
   fprintf(stderr, "\n fatal error in Chv_setRealEntry(%p,%d,%d,%e)"
           "\n bad type %d, not SPOOLES_REAL\n", 
           chv, irow, jcol, value, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_setRealEntry(%p,%d,%d,%e)"
           "\n bad symflag %d"
           "\n must be SPOOLES_SYMMETRIC of SPOOLES_NONSYMMETRIC\n",
           chv, irow, jcol, value, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, "\n fatal error in Chv_setRealEntry(%p,%d,%d,%e)"
           "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
           chv, irow, jcol, value, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow < nD || jcol < nD ) {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( CHV_IS_SYMMETRIC(chv) && off < 0 ) {
      off = -off ;
   }
   base = Chv_diagLocation(chv, ichv) ;
   base[off] = value ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   fill (*pReal,*pImag) with entry (irow, jcol)
 
   created -- 98apr30, cca
   --------------------------------------------
*/
void
Chv_setComplexEntry (
   Chv      *chv,
   int      irow,
   int      jcol,
   double   real,
   double   imag
) {
int      ichv, ncol, nD, nL, nrow, nU, off ;
double   *base ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || irow < 0 || jcol < 0 ) {
   fprintf(stderr, 
           "\n fatal error in Chv_setComplexEntry(%p,%d,%d,%e,%e)"
           "\n bad input\n", chv, irow, jcol, real, imag) ;
   exit(-1) ;
}
if ( ! CHV_IS_COMPLEX(chv) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_setComplexEntry(%p,%d,%d,%e,%e)"
           "\n bad type %d, not SPOOLES_COMPLEX\n", 
           chv, irow, jcol, real, imag, chv->type) ;
   exit(-1) ;
}
if ( ! (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) 
        || CHV_IS_NONSYMMETRIC(chv)) ) {
   fprintf(stderr, 
           "\n fatal error in Chv_setComplexEntry(%p,%d,%d,%e,%e)"
           "\n bad symflag %d"
           "\n not SPOOLES_SYMMETRIC, SPOOLES_HERMITIAN"
           "\n or SPOOLES_NONSYMMETRIC \n",
           chv, irow, jcol, real, imag, chv->symflag) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
ncol = nD + nU ;
if ( CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv) ) {
   nrow = ncol ;
} else {
   nrow = nD + nL ;
}
if ( irow >= nrow || jcol >= ncol ) {
   fprintf(stderr, 
           "\n fatal error in Chv_setComplexEntry(%p,%d,%d,%e,%e)"
           "\n irow = %d, jcol = %d, nrow = %d, ncol = %d\n",
           chv, irow, jcol, real, imag, irow, jcol, nrow, ncol) ;
   exit(-1) ;
}
if ( irow < nD || jcol < nD ) {
   ichv = (irow <= jcol) ? irow : jcol ;
   off  = jcol - irow ;
   if ( off < 0 && (CHV_IS_SYMMETRIC(chv) || CHV_IS_HERMITIAN(chv)) ) {
      off = -off ;
   }
   base   = Chv_diagLocation(chv, ichv) ;
   base[2*off]   = real ;
   base[2*off+1] = imag ;
}
return ; }

/*--------------------------------------------------------------------*/
