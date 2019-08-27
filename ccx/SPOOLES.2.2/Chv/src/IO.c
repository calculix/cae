/*  IO.c  */

#include "../Chv.h"

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------
   purpose -- to write the object to a file
              in human readable form

   created -- 98apr30, cca
   ----------------------------------------
*/
void
Chv_writeForHumanEye (
   Chv    *chv,
   FILE   *fp
) {
A2    mtx ;
int   ierr, ncol, nD, nL, nrow, nU ;
int   *colind, *rowind ; 
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_writeForHumanEye(%p,%p)"
           "\n bad input\n", chv, fp) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
fprintf(fp, 
       "\n Chv object at address %p"
       "\n id = %d, nD = %d, nL = %d, nU = %d, type = %d, symflag = %d",
       chv, chv->id, nD, nL, nU, chv->type, chv->symflag) ;
if ( CHV_IS_REAL(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
      fprintf(fp, "\n chv is real and symmetric") ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      fprintf(fp, "\n chv is real and nonsymmetric") ;
   } else {
      fprintf(fp, "\n chv has unknown symmetry type %d", chv->symflag) ;
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   if ( CHV_IS_SYMMETRIC(chv) ) {
      fprintf(fp, "\n chv is complex and symmetric") ;
   } else if ( CHV_IS_HERMITIAN(chv) ) {
      fprintf(fp, "\n chv is complex and hermitian") ;
   } else if ( CHV_IS_NONSYMMETRIC(chv) ) {
      fprintf(fp, "\n chv is complex and nonsymmetric") ;
   } else {
      fprintf(fp, "\n chv has unknown symmetry type %d", chv->symflag) ;
   }
} else {
   fprintf(fp, "\n chv has unknown type %d", chv->type) ;
}
Chv_rowIndices(chv, &nrow, &rowind) ;
if ( nrow > 0 && rowind != NULL ) {
   fprintf(fp, "\n chv's row indices at %p", rowind) ;
   IVfp80(fp, nrow, rowind, 80, &ierr) ;
}
Chv_columnIndices(chv, &ncol, &colind) ;
if ( ncol > 0 && colind != NULL ) {
   fprintf(fp, "\n chv's column indices at %p", colind) ;
   IVfp80(fp, ncol, colind, 80, &ierr) ;
}
/*
   --------------------
   load the (1,1) block
   --------------------
*/
A2_setDefaultFields(&mtx) ;
Chv_fill11block(chv, &mtx) ;
fprintf(fp, "\n (1,1) block") ;
A2_writeForHumanEye(&mtx, fp) ;
if ( nU > 0 ) {
/*
   --------------------
   load the (1,2) block
   --------------------
*/
   Chv_fill12block(chv, &mtx) ;
   fprintf(fp, "\n (1,2) block") ;
   A2_writeForHumanEye(&mtx, fp) ;
}
if ( nL > 0 && CHV_IS_NONSYMMETRIC(chv) == 1 ) {
/*
   --------------------
   load the (2,1) block
   --------------------
*/
   Chv_fill21block(chv, &mtx) ;
   fprintf(fp, "\n (2,1) block") ;
   A2_writeForHumanEye(&mtx, fp) ;
}
A2_clearData(&mtx) ;

return ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   purpose -- write out the entries in matlab style

   created -- 98apr30, cca
   ------------------------------------------------
*/
void
Chv_writeForMatlab (
   Chv    *chv,
   char   *chvname,
   FILE   *fp
) {
int      irow, jcol, ncol, nD, nL, nrow, nU ;
int      *colind, *rowind ;
/*
   ---------------
   check the input
   ---------------
*/
if ( chv == NULL || chvname == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in Chv_writeForMatlab(%p,%p,%p)"
           "\n bad input\n", chv, chvname, fp) ;
   exit(-1) ;
}
if ( ! (CHV_IS_REAL(chv) || CHV_IS_COMPLEX(chv)) ) {
   fprintf(stderr, "\n fatal error in Chv_writeForMatlab(%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           chv, chvname, fp, chv->type) ;
   exit(-1) ;
}
Chv_dimensions(chv, &nD, &nL, &nU) ;
Chv_rowIndices(chv, &nrow, &rowind) ;
Chv_columnIndices(chv, &ncol, &colind) ;
if ( CHV_IS_REAL(chv) ) {
   double   value ;
/*
   -------------------------
   write out the (1,1) block
   -------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      for ( jcol = 0 ; jcol < nD ; jcol++ ) {
         Chv_realEntry(chv, irow, jcol, &value) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e ;",
                 chvname, 1+rowind[irow], 1+colind[jcol], value) ;
      }
   }
/*
   -------------------------
   write out the (1,2) block
   -------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      for ( jcol = nD ; jcol < ncol ; jcol++ ) {
         Chv_realEntry(chv, irow, jcol, &value) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e ;",
                 chvname, 1+rowind[irow], 1+colind[jcol], value) ;
      }
   }
/*
   -------------------------
   write out the (2,1) block
   -------------------------
*/
   for ( irow = nD ; irow < nrow ; irow++ ) {
      for ( jcol = 0 ; jcol < nD ; jcol++ ) {
         Chv_realEntry(chv, irow, jcol, &value) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e ;",
                 chvname, 1+rowind[irow], 1+colind[jcol], value) ;
      }
   }
} else if ( CHV_IS_COMPLEX(chv) ) {
   double   imag, real ;
/*
   -------------------------
   write out the (1,1) block
   -------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      for ( jcol = 0 ; jcol < nD ; jcol++ ) {
         Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;",
                 chvname, 1+rowind[irow], 1+colind[jcol],
                 real, imag) ;
      }
   }
/*
   -------------------------
   write out the (1,2) block
   -------------------------
*/
   for ( irow = 0 ; irow < nD ; irow++ ) {
      for ( jcol = nD ; jcol < ncol ; jcol++ ) {
         Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;",
                 chvname, 1+rowind[irow], 1+colind[jcol],
                 real, imag) ;
      }
   }
/*
   -------------------------
   write out the (2,1) block
   -------------------------
*/
   for ( irow = nD ; irow < nrow ; irow++ ) {
      for ( jcol = 0 ; jcol < nD ; jcol++ ) {
         Chv_complexEntry(chv, irow, jcol, &real, &imag) ;
         fprintf(fp, "\n %s(%d,%d) = %20.12e + %20.12e*i;",
                 chvname, 1+rowind[irow], 1+colind[jcol],
                 real, imag) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
