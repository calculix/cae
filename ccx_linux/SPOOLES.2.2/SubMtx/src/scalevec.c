/*  scalevec.c  */

#include "../SubMtx.h"

/*--------------------------------------------------------------------*/
static void diagonal_scale3vec ( SubMtx *mtxA, double y0[], double y1[],
   double y2[], double x0[], double x1[], double x2[] ) ;
static void diagonal_scale2vec ( SubMtx *mtxA, double y0[], double y1[],
   double x0[], double x1[] ) ;
static void diagonal_scale1vec ( SubMtx *mtxA, 
   double y0[], double x0[] ) ;
static void block_diagonal_sym_scale3vec ( SubMtx *mtxA, double y0[],
   double y1[], double y2[], double x0[], double x1[], double x2[] ) ;
static void block_diagonal_sym_scale2vec ( SubMtx *mtxA, double y0[],
   double y1[], double x0[], double x1[] ) ;
static void block_diagonal_sym_scale1vec ( SubMtx *mtxA, double y0[],
   double x0[] ) ;
static void block_diagonal_herm_scale3vec ( SubMtx *mtxA, double y0[],
   double y1[], double y2[], double x0[], double x1[], double x2[] ) ;
static void block_diagonal_herm_scale2vec ( SubMtx *mtxA, double y0[],
   double y1[], double x0[], double x1[] ) ;
static void block_diagonal_herm_scale1vec ( SubMtx *mtxA, double y0[],
   double x0[] ) ;
/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- compute [y0] := A * [x0]

   created -- 98apr17, cca
   -----------------------------------
*/
void
SubMtx_scale1vec (
   SubMtx   *mtxA,
   double   y0[],
   double   x0[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if ( mtxA == NULL || y0 == NULL || x0 == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_scale1vec(%p,%p,%p)"
           "\n bad input\n", mtxA, y0, x0) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtxA) || SUBMTX_IS_COMPLEX(mtxA)) ) {
   fprintf(stderr, "\n fatal error in SubMtx_scale1vec(%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtxA, y0, x0, mtxA->type) ;
   exit(-1) ;
}
switch ( mtxA->mode ) {
case SUBMTX_DIAGONAL :
   diagonal_scale1vec(mtxA, y0, x0) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
   block_diagonal_sym_scale1vec(mtxA, y0, x0) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_HERM :
   if ( ! SUBMTX_IS_COMPLEX(mtxA) ) {
      fprintf(stderr, "\n fatal error in SubMtx_scale1vec(%p,%p,%p)"
              "\n hermitian matrix, type %d is not SPOOLES_COMPLEX\n", 
              mtxA, y0, x0, mtxA->type) ;
      exit(-1) ;
   }
   block_diagonal_herm_scale1vec(mtxA, y0, x0) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_scale1vec()"
           "\n matrix mode not supported"
           "\n must be SUBMTX_DIAGONAL,"
           "\n      or SUBMTX_BLOCK_DIAGONAL_SYM"
           "\n      or SUBMTX_BLOCK_DIAGONAL_HERM\n") ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -------------------------------------
   purpose -- compute [y0 y1] := [x0 x1]

   created -- 98apr17, cca
   -------------------------------------
*/
void
SubMtx_scale2vec (
   SubMtx     *mtxA,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxA == NULL || y0 == NULL || y1 == NULL 
   || x0 == NULL || x1 == NULL ) {
   fprintf(stderr, "\n fatal error in SubMtx_scale2vec(%p,%p,%p,%p,%p)"
           "\n bad input\n", mtxA, y0, y1, x0, x1) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtxA) || SUBMTX_IS_COMPLEX(mtxA)) ) {
   fprintf(stderr, "\n fatal error in SubMtx_scale2vec(%p,%p,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtxA, y0, y1, x0, x1, mtxA->type) ;
   exit(-1) ;
}
switch ( mtxA->mode ) {
case SUBMTX_DIAGONAL :
   diagonal_scale2vec(mtxA, y0, y1, x0, x1) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
   block_diagonal_sym_scale2vec(mtxA, y0, y1, x0, x1) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_HERM :
   if ( ! SUBMTX_IS_COMPLEX(mtxA) ) {
      fprintf(stderr, 
              "\n fatal error in SubMtx_scale2vec(%p,%p,%p,%p,%p)"
              "\n hermitian matrix, type %d is not SPOOLES_COMPLEX\n", 
              mtxA, y0, y1, x0, x1, mtxA->type) ;
      exit(-1) ;
   }
   block_diagonal_herm_scale2vec(mtxA, y0, y1, x0, x1) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_scale2vec()"
           "\n matrix type not supported"
           "\n must be SUBMTX_DIAGONAL,"
           "\n      or SUBMTX_BLOCK_DIAGONAL_SYM"
           "\n      or SUBMTX_BLOCK_DIAGONAL_HERM\n") ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- compute [y0 y1 y2] := A * [x0 x1 x2]

   created -- 98apr17, cca
   -----------------------------------------------
*/
void
SubMtx_scale3vec (
   SubMtx     *mtxA,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxA == NULL || y0 == NULL || y1 == NULL || y2 == NULL
   || x0 == NULL || x1 == NULL || x2 == NULL ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_scale3vec(%p,%p,%p,%p,%p,%p,%p)"
           "\n bad input\n", mtxA, y0, y1, y2, x0, x1, x2) ;
   exit(-1) ;
}
if ( ! (SUBMTX_IS_REAL(mtxA) || SUBMTX_IS_COMPLEX(mtxA)) ) {
   fprintf(stderr, 
           "\n fatal error in SubMtx_scale3vec(%p,%p,%p,%p,%p,%p,%p)"
           "\n bad type %d, must be SPOOLES_REAL or SPOOLES_COMPLEX\n", 
           mtxA, y0, y1, y2, x0, x1, x2, mtxA->type) ;
   exit(-1) ;
}
switch ( mtxA->mode ) {
case SUBMTX_DIAGONAL :
   diagonal_scale3vec(mtxA, y0, y1, y2, x0, x1, x2) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_SYM :
   block_diagonal_sym_scale3vec(mtxA, y0, y1, y2, x0, x1, x2) ;
   break ;
case SUBMTX_BLOCK_DIAGONAL_HERM :
   if ( ! SUBMTX_IS_COMPLEX(mtxA) ) {
      fprintf(stderr, 
              "\n fatal error in SubMtx_scale3vec(%p,%p,%p,%p,%p,%p,%p)"
              "\n hermitian matrix, type %d is not SPOOLES_COMPLEX\n", 
              mtxA, y0, y1, y2, x0, x1, x2, mtxA->type) ;
      exit(-1) ;
   }
   block_diagonal_herm_scale3vec(mtxA, y0, y1, y2, x0, x1, x2) ;
   break ;
default :
   fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
           "\n matrix type not supported"
           "\n must be SUBMTX_DIAGONAL,"
           "\n      or SUBMTX_BLOCK_DIAGONAL_SYM"
           "\n      or SUBMTX_BLOCK_DIAGONAL_HERM\n") ;
   exit(-1) ;
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   purpose -- compute [y0 y1 y2] := A * [x0 x1 x2]
              where A is diagonal

   created -- 98apr17, cca
   -----------------------------------------------
*/
static void
diagonal_scale3vec (
   SubMtx   *mtxA,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
double   *entriesA ;
int      nrowA ;

SubMtx_diagonalInfo(mtxA, &nrowA, &entriesA) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   double   a ;
   int      irowA ;
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      a  = entriesA[irowA] ;
      y0[irowA] = a*x0[irowA] ; 
      y1[irowA] = a*x1[irowA] ; 
      y2[irowA] = a*x2[irowA] ; 
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   double   ai, ar, xi0, xi1, xi2, xr0, xr1, xr2 ;
   int      iloc, irowA, rloc ;
   for ( irowA = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         irowA++, rloc += 2, iloc += 2 ) {
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      xr1 = x1[rloc] ; xi1 = x1[iloc] ;
      xr2 = x2[rloc] ; xi2 = x2[iloc] ;
      ar  = entriesA[rloc] ; ai = entriesA[iloc] ;
      y0[rloc] = ar*xr0 - ai*xi0 ; y0[iloc] = ar*xi0 + ai*xr0 ;
      y1[rloc] = ar*xr1 - ai*xi1 ; y1[iloc] = ar*xi1 + ai*xr1 ;
      y2[rloc] = ar*xr2 - ai*xi2 ; y2[iloc] = ar*xi2 + ai*xr2 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute [y0 y1] := A * [x0 x1]
              where A is diagonal

   created -- 98apr17, cca
   -----------------------------------------
*/
static void
diagonal_scale2vec (
   SubMtx   *mtxA,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[]
) {
double   *entriesA ;
int      nrowA ;

SubMtx_diagonalInfo(mtxA, &nrowA, &entriesA) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   double   a ;
   int      irowA ;
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      a  = entriesA[irowA] ;
      y0[irowA] = a*x0[irowA] ; 
      y1[irowA] = a*x1[irowA] ; 
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   double   ai, ar, xi0, xi1, xr0, xr1 ;
   int      iloc, irowA, rloc ;
   for ( irowA = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         irowA++, rloc += 2, iloc += 2 ) {
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      xr1 = x1[rloc] ; xi1 = x1[iloc] ;
      ar  = entriesA[rloc] ; ai = entriesA[iloc] ;
      y0[rloc] = ar*xr0 - ai*xi0 ; y0[iloc] = ar*xi0 + ai*xr0 ;
      y1[rloc] = ar*xr1 - ai*xi1 ; y1[iloc] = ar*xi1 + ai*xr1 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   purpose -- compute [y0] := A * [x0]
              where A is diagonal

   created -- 98apr17, cca
   -----------------------------------
*/
static void
diagonal_scale1vec (
   SubMtx   *mtxA,
   double   y0[],
   double   x0[]
) {
double   *entriesA ;
int      nrowA ;

SubMtx_diagonalInfo(mtxA, &nrowA, &entriesA) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   double   a ;
   int      irowA ;
   for ( irowA = 0 ; irowA < nrowA ; irowA++ ) {
      a  = entriesA[irowA] ;
      y0[irowA] = a*x0[irowA] ; 
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   double   ai, ar, xi0, xr0 ;
   int      iloc, irowA, rloc ;
   for ( irowA = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         irowA++, rloc += 2, iloc += 2 ) {
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      ar  = entriesA[rloc] ; ai = entriesA[iloc] ;
      y0[rloc] = ar*xr0 - ai*xi0 ; y0[iloc] = ar*xi0 + ai*xr0 ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0 y1 y2] := A * [x0 x1 x2]
              where A is block diagonal and symmetric

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_sym_scale3vec (
   SubMtx   *mtxA,
   double   y0[],
   double   y1[],
   double   y2[],
   double   x0[],
   double   x1[],
   double   x2[]
) {
double   *entries ;
int      nentA, nrowA ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   int      ipivot, irowA, kk ;

   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   a00, x00, x01, x02 ;
   
         a00 = entries[kk] ; 
         x00 = x0[irowA] ; 
         x01 = x1[irowA] ; 
         x02 = x2[irowA] ; 
         y0[irowA] = a00*x00 ;
         y1[irowA] = a00*x01 ;
         y2[irowA] = a00*x02 ;
         kk++, irowA++ ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   a00, a01, a11, 
                  x00, x01, x02, x10, x11, x12 ;
   
         a00 = entries[kk]   ; 
         a01 = entries[kk+1] ; 
         a11 = entries[kk+2] ; 
         x00 = x0[irowA]     ;
         x01 = x1[irowA]     ;
         x02 = x2[irowA]     ;
         x10 = x0[irowA+1]   ;
         x11 = x1[irowA+1]   ;
         x12 = x2[irowA+1]   ;
         y0[irowA]   = a00*x00 + a01*x10 ;
         y1[irowA]   = a00*x01 + a01*x11 ;
         y2[irowA]   = a00*x02 + a01*x12 ;
         y0[irowA+1] = a01*x00 + a11*x10 ;
         y1[irowA+1] = a01*x01 + a11*x11 ;
         y2[irowA+1] = a01*x02 + a11*x12 ;
         kk += 3, irowA += 2 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   int      iloc, ipivot, irowA, kk, rloc ;

   for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   ai00, ar00, xi0, xi1, xi2, xr0, xr1, xr2 ;
   
         ar00 = entries[kk] ; ai00 = entries[kk+1] ;
         xr0 = x0[rloc] ; xi0 = x0[iloc] ;
         xr1 = x1[rloc] ; xi1 = x1[iloc] ;
         xr2 = x2[rloc] ; xi2 = x2[iloc] ;
         y0[rloc] = ar00*xr0 - ai00*xi0; y0[iloc] = ar00*xi0 + ai00*xr0;
         y1[rloc] = ar00*xr1 - ai00*xi1; y1[iloc] = ar00*xi1 + ai00*xr1;
         y2[rloc] = ar00*xr2 - ai00*xi2; y2[iloc] = ar00*xi2 + ai00*xr2;
         kk += 2, irowA++, rloc += 2, iloc += 2 ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   ai00, ai01, ai11, ar00, ar01, ar11,
                  xi00, xi01, xi02, xi10, xi11, xi12,
                  xr00, xr01, xr02, xr10, xr11, xr12 ;
         int      iloc0, iloc1, rloc0, rloc1 ;
   
         ar00 = entries[kk]   ; ai00 = entries[kk+1] ;
         ar01 = entries[kk+2] ; ai01 = entries[kk+3] ;
         ar11 = entries[kk+4] ; ai11 = entries[kk+5] ;
         rloc0 = rloc     ; iloc0 = iloc ;
         rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
         xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
         xr01 = x1[rloc0] ; xi01 = x1[iloc0] ;
         xr02 = x2[rloc0] ; xi02 = x2[iloc0] ;
         xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
         xr11 = x1[rloc1] ; xi11 = x1[iloc1] ;
         xr12 = x2[rloc1] ; xi12 = x2[iloc1] ;
         y0[rloc0] = ar00*xr00 - ai00*xi00 + ar01*xr10 - ai01*xi10 ;
         y0[iloc0] = ar00*xi00 + ai00*xr00 + ar01*xi10 + ai01*xr10 ;
         y1[rloc0] = ar00*xr01 - ai00*xi01 + ar01*xr11 - ai01*xi11 ;
         y1[iloc0] = ar00*xi01 + ai00*xr01 + ar01*xi11 + ai01*xr11 ;
         y2[rloc0] = ar00*xr02 - ai00*xi02 + ar01*xr12 - ai01*xi12 ;
         y2[iloc0] = ar00*xi02 + ai00*xr02 + ar01*xi12 + ai01*xr12 ;
         y0[rloc1] = ar01*xr00 - ai01*xi00 + ar11*xr10 - ai11*xi10 ;
         y0[iloc1] = ar01*xi00 + ai01*xr00 + ar11*xi10 + ai11*xr10 ;
         y1[rloc1] = ar01*xr01 - ai01*xi01 + ar11*xr11 - ai11*xi11 ;
         y1[iloc1] = ar01*xi01 + ai01*xr01 + ar11*xi11 + ai11*xr11 ;
         y2[rloc1] = ar01*xr02 - ai01*xi02 + ar11*xr12 - ai11*xi12 ;
         y2[iloc1] = ar01*xi02 + ai01*xr02 + ar11*xi12 + ai11*xr12 ;
         kk += 6, irowA += 2 ; rloc += 4 ; iloc += 4 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0 y1] := A * [x0 x1]
              where A is block diagonal and symmetric

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_sym_scale2vec (
   SubMtx   *mtxA,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[]
) {
double   *entries ;
int      nentA, nrowA ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   int      ipivot, irowA, kk ;

   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   a00, x00, x01 ;
   
         a00 = entries[kk] ; 
         x00 = x0[irowA] ; 
         x01 = x1[irowA] ; 
         y0[irowA] = a00*x00 ;
         y1[irowA] = a00*x01 ;
         kk++, irowA++ ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   a00, a01, a11, x00, x01, x10, x11 ;
   
         a00 = entries[kk]   ; 
         a01 = entries[kk+1] ; 
         a11 = entries[kk+2] ; 
         x00 = x0[irowA]     ;
         x01 = x1[irowA]     ;
         x10 = x0[irowA+1]   ;
         x11 = x1[irowA+1]   ;
         y0[irowA]   = a00*x00 + a01*x10 ;
         y1[irowA]   = a00*x01 + a01*x11 ;
         y0[irowA+1] = a01*x00 + a11*x10 ;
         y1[irowA+1] = a01*x01 + a11*x11 ;
         kk += 3, irowA += 2 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   int      iloc, ipivot, irowA, kk, rloc ;

   for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   ai00, ar00, xi0, xi1, xr0, xr1 ;
   
         ar00 = entries[kk] ; ai00 = entries[kk+1] ;
         xr0 = x0[rloc] ; xi0 = x0[iloc] ;
         xr1 = x1[rloc] ; xi1 = x1[iloc] ;
         y0[rloc] = ar00*xr0 - ai00*xi0; y0[iloc] = ar00*xi0 + ai00*xr0;
         y1[rloc] = ar00*xr1 - ai00*xi1; y1[iloc] = ar00*xi1 + ai00*xr1;
         kk += 2, irowA++, rloc += 2, iloc += 2 ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   ai00, ai01, ai11, ar00, ar01, ar11,
                  xi00, xi01, xi10, xi11, xr00, xr01, xr10, xr11 ;
         int      iloc0, iloc1, rloc0, rloc1 ;
   
         ar00 = entries[kk]   ; ai00 = entries[kk+1] ;
         ar01 = entries[kk+2] ; ai01 = entries[kk+3] ;
         ar11 = entries[kk+4] ; ai11 = entries[kk+5] ;
         rloc0 = rloc     ; iloc0 = iloc ;
         rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
         xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
         xr01 = x1[rloc0] ; xi01 = x1[iloc0] ;
         xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
         xr11 = x1[rloc1] ; xi11 = x1[iloc1] ;
         y0[rloc0] = ar00*xr00 - ai00*xi00 + ar01*xr10 - ai01*xi10 ;
         y0[iloc0] = ar00*xi00 + ai00*xr00 + ar01*xi10 + ai01*xr10 ;
         y1[rloc0] = ar00*xr01 - ai00*xi01 + ar01*xr11 - ai01*xi11 ;
         y1[iloc0] = ar00*xi01 + ai00*xr01 + ar01*xi11 + ai01*xr11 ;
         y0[rloc1] = ar01*xr00 - ai01*xi00 + ar11*xr10 - ai11*xi10 ;
         y0[iloc1] = ar01*xi00 + ai01*xr00 + ar11*xi10 + ai11*xr10 ;
         y1[rloc1] = ar01*xr01 - ai01*xi01 + ar11*xr11 - ai11*xi11 ;
         y1[iloc1] = ar01*xi01 + ai01*xr01 + ar11*xi11 + ai11*xr11 ;
         kk += 6, irowA += 2 ; rloc += 4 ; iloc += 4 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale2vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0] := A * [x0]
              where A is block diagonal and symmetric

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_sym_scale1vec (
   SubMtx   *mtxA,
   double   y0[],
   double   x0[]
) {
double   *entries ;
int      nentA, nrowA ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
if ( SUBMTX_IS_REAL(mtxA) ) {
   int      ipivot, irowA, kk ;

   for ( irowA = ipivot = kk = 0 ; irowA < nrowA ; ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   a00, x00 ;
   
         a00 = entries[kk] ; 
         x00 = x0[irowA] ; 
         y0[irowA] = a00*x00 ;
         kk++, irowA++ ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   a00, a01, a11, x00, x10 ;
   
         a00 = entries[kk]   ; 
         a01 = entries[kk+1] ; 
         a11 = entries[kk+2] ; 
         x00 = x0[irowA]     ;
         x10 = x0[irowA+1]   ;
         y0[irowA]   = a00*x00 + a01*x10 ;
         y0[irowA+1] = a01*x00 + a11*x10 ;
         kk += 3, irowA += 2 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
} else if ( SUBMTX_IS_COMPLEX(mtxA) ) {
   int      iloc, ipivot, irowA, kk, rloc ;

   for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
         irowA < nrowA ; 
         ipivot++ ) {
      if ( pivotsizes[ipivot] == 1 ) {
         double   ai00, ar00, xi0, xr0 ;
   
         ar00 = entries[kk] ; ai00 = entries[kk+1] ;
         xr0 = x0[rloc] ; xi0 = x0[iloc] ;
         y0[rloc] = ar00*xr0 - ai00*xi0; y0[iloc] = ar00*xi0 + ai00*xr0;
         kk += 2, irowA++, rloc += 2, iloc += 2 ;
      } else if ( pivotsizes[ipivot] == 2 ) {
         double   ai00, ai01, ai11, ar00, ar01, ar11,
                  xi00, xi10, xr00, xr10 ;
         int      iloc0, iloc1, rloc0, rloc1 ;
   
         ar00 = entries[kk]   ; ai00 = entries[kk+1] ;
         ar01 = entries[kk+2] ; ai01 = entries[kk+3] ;
         ar11 = entries[kk+4] ; ai11 = entries[kk+5] ;
         rloc0 = rloc     ; iloc0 = iloc ;
         rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
         xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
         xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
         y0[rloc0] = ar00*xr00 - ai00*xi00 + ar01*xr10 - ai01*xi10 ;
         y0[iloc0] = ar00*xi00 + ai00*xr00 + ar01*xi10 + ai01*xr10 ;
         y0[rloc1] = ar01*xr00 - ai01*xi00 + ar11*xr10 - ai11*xi10 ;
         y0[iloc1] = ar01*xi00 + ai01*xr00 + ar11*xi10 + ai11*xr10 ;
         kk += 6, irowA += 2 ; rloc += 4 ; iloc += 4 ;
      } else {
         fprintf(stderr, "\n fatal error in SubMtx_scale1vec()"
                 "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
         exit(-1) ;
      }
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0 y1 y2] := A * [x0 x1 x2]
              where A is block diagonal and hermitian

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_herm_scale3vec (
   SubMtx     *mtxA,
   double     y0[],
   double     y1[],
   double     y2[],
   double     x0[],
   double     x1[],
   double     x2[]
) {
double   *entries ;
int      iloc, ipivot, irowA, kk, nentA, nrowA, rloc ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
      irowA < nrowA ; 
      ipivot++ ) {
   if ( pivotsizes[ipivot] == 1 ) {
      double   ar00, ai00, xi0, xi1, xi2, xr0, xr1, xr2 ;

/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ; 
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      xr1 = x1[rloc] ; xi1 = x1[iloc] ;
      xr2 = x2[rloc] ; xi2 = x2[iloc] ;
      y0[rloc] = ar00*xr0 ; y0[iloc] = ar00*xi0 ;
      y1[rloc] = ar00*xr1 ; y1[iloc] = ar00*xi1 ;
      y2[rloc] = ar00*xr2 ; y2[iloc] = ar00*xi2 ;
      irowA++, rloc += 2, iloc += 2 ;
   } else if ( pivotsizes[ipivot] == 2 ) {
      double   ai00, ai01, ai11, ar00, ar01, ar11,
               xi00, xi01, xi02, xi10, xi11, xi12,
               xr00, xr01, xr02, xr10, xr11, xr12 ;
      int      iloc0, iloc1, rloc0, rloc1 ;

      rloc0 = rloc     ; iloc0 = iloc     ;
      rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = entries[kk++] ;
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = 0.0 ; kk++ ;
      xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
      xr01 = x1[rloc0] ; xi01 = x1[iloc0] ;
      xr02 = x2[rloc0] ; xi02 = x2[iloc0] ;
      xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
      xr11 = x1[rloc1] ; xi11 = x1[iloc1] ;
      xr12 = x2[rloc1] ; xi12 = x2[iloc1] ;
      y0[rloc0] = ar00*xr00 + ar01*xr10 - ai01*xi10 ;
      y0[iloc0] = ar00*xi00 + ar01*xi10 + ai01*xr10 ;
      y1[rloc0] = ar00*xr01 + ar01*xr11 - ai01*xi11 ;
      y1[iloc0] = ar00*xi01 + ar01*xi11 + ai01*xr11 ;
      y2[rloc0] = ar00*xr02 + ar01*xr12 - ai01*xi12 ;
      y2[iloc0] = ar00*xi02 + ar01*xi12 + ai01*xr12 ;
      y0[rloc1] = ar01*xr00 + ai01*xi00 + ar11*xr10 ;
      y0[iloc1] = ar01*xi00 - ai01*xr00 + ar11*xi10 ;
      y1[rloc1] = ar01*xr01 + ai01*xi01 + ar11*xr11 ;
      y1[iloc1] = ar01*xi01 - ai01*xr01 + ar11*xi11 ;
      y2[rloc1] = ar01*xr02 + ai01*xi02 + ar11*xr12 ;
      y2[iloc1] = ar01*xi02 - ai01*xr02 + ar11*xi12 ;
      irowA += 2 ; rloc += 4 ; iloc += 4 ;
   } else {
      fprintf(stderr, "\n fatal error in SubMtx_scale3vec()"
              "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
      exit(-1) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0 y1] := A * [x0 x1]
              where A is block diagonal and hermitian

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_herm_scale2vec (
   SubMtx   *mtxA,
   double   y0[],
   double   y1[],
   double   x0[],
   double   x1[]
) {
double   *entries ;
int      iloc, ipivot, irowA, kk, nentA, nrowA, rloc ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
      irowA < nrowA ; 
      ipivot++ ) {
   if ( pivotsizes[ipivot] == 1 ) {
      double   ai00, ar00, xi0, xi1, xr0, xr1 ;

/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ;
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      xr1 = x1[rloc] ; xi1 = x1[iloc] ;
      y0[rloc] = ar00*xr0 - ai00*xi0 ; y0[iloc] = ar00*xi0 + ai00*xr0 ;
      y1[rloc] = ar00*xr1 - ai00*xi1 ; y1[iloc] = ar00*xi1 + ai00*xr1 ;
      irowA++, rloc += 2, iloc += 2 ;
   } else if ( pivotsizes[ipivot] == 2 ) {
      double   ai00, ai01, ai11, ar00, ar01, ar11,
               xi00, xi01, xi10, xi11, xr00, xr01, xr10, xr11 ;
      int      iloc0, iloc1, rloc0, rloc1 ;

      rloc0 = rloc     ; iloc0 = iloc     ;
      rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = entries[kk++] ;
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = 0.0 ; kk++ ;
      xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
      xr01 = x1[rloc0] ; xi01 = x1[iloc0] ;
      xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
      xr11 = x1[rloc1] ; xi11 = x1[iloc1] ;
      y0[rloc0] = ar00*xr00 + ar01*xr10 - ai01*xi10 ;
      y0[iloc0] = ar00*xi00 + ar01*xi10 + ai01*xr10 ;
      y1[rloc0] = ar00*xr01 + ar01*xr11 - ai01*xi11 ;
      y1[iloc0] = ar00*xi01 + ar01*xi11 + ai01*xr11 ;
      y0[rloc1] = ar01*xr00 + ai01*xi00 + ar11*xr10 ;
      y0[iloc1] = ar01*xi00 - ai01*xr00 + ar11*xi10 ;
      y1[rloc1] = ar01*xr01 + ai01*xi01 + ar11*xr11 ;
      y1[iloc1] = ar01*xi01 - ai01*xr01 + ar11*xi11 ;
      irowA += 2 ; rloc += 4 ; iloc += 4 ;
   } else {
      fprintf(stderr, "\n fatal error in SubMtx_scale2vec()"
              "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
      exit(-1) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   purpose -- compute [y0] := A * [x0]
              where A is block diagonal and hermitian

   created -- 98apr17, cca
   --------------------------------------------------
*/
static void
block_diagonal_herm_scale1vec (
   SubMtx   *mtxA,
   double   y0[],
   double   x0[]
) {
double   *entries ;
int      iloc, ipivot, irowA, kk, nentA, nrowA, rloc ;
int      *pivotsizes ;

SubMtx_blockDiagonalInfo(mtxA, &nrowA, &nentA, &pivotsizes, &entries) ;
for ( irowA = ipivot = kk = rloc = 0, iloc = 1 ; 
      irowA < nrowA ; 
      ipivot++ ) {
   if ( pivotsizes[ipivot] == 1 ) {
      double   ai00, ar00, xi0, xr0 ;

/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ;
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      xr0 = x0[rloc] ; xi0 = x0[iloc] ;
      y0[rloc] = ar00*xr0 - ai00*xi0 ; y0[iloc] = ar00*xi0 + ai00*xr0 ;
      irowA++, rloc += 2, iloc += 2 ;
   } else if ( pivotsizes[ipivot] == 2 ) {
      double   ai00, ai01, ai11, ar00, ar01, ar11,
               xi00, xi10, xr00, xr10 ;
      int      iloc0, iloc1, rloc0, rloc1 ;

      rloc0 = rloc     ; iloc0 = iloc     ;
      rloc1 = rloc + 2 ; iloc1 = iloc + 2 ;
/*
      ar00 = entries[kk++] ; ai00 = entries[kk++] ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = entries[kk++] ;
*/
      ar00 = entries[kk++] ; ai00 = 0.0 ; kk++ ;
      ar01 = entries[kk++] ; ai01 = entries[kk++] ;
      ar11 = entries[kk++] ; ai11 = 0.0 ; kk++ ;
      xr00 = x0[rloc0] ; xi00 = x0[iloc0] ;
      xr10 = x0[rloc1] ; xi10 = x0[iloc1] ;
      y0[rloc0] = ar00*xr00 + ar01*xr10 - ai01*xi10 ;
      y0[iloc0] = ar00*xi00 + ar01*xi10 + ai01*xr10 ;
      y0[rloc1] = ar01*xr00 + ai01*xi00 + ar11*xr10 ;
      y0[iloc1] = ar01*xi00 - ai01*xr00 + ar11*xi10 ;
      irowA += 2 ; rloc += 4 ; iloc += 4 ;
   } else {
      fprintf(stderr, "\n fatal error in SubMtx_scale1vec()"
              "\n pivotsizes[%d] = %d", ipivot, pivotsizes[ipivot]) ;
      exit(-1) ;
   }
}
return ; }

/*--------------------------------------------------------------------*/
