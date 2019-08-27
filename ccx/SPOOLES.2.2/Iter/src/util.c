#include "../Iter.h"

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   return the frobenius norm of a Dense matrix

   created -- 98nov11, dkw
   -------------------------------------------
*/
double
DenseMtx_frobNorm (
   DenseMtx *mtx
) {
double fnorm;
A2 a2;
A2_setDefaultFields(&a2) ;
DenseMtx_setA2(mtx, &a2);
fnorm = A2_frobNorm(&a2) ;
return (fnorm); }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   return the two-norm of column jcol from a Dense matrix

   created -- 98nov11, dkw
   ------------------------------------------------------
*/
double
DenseMtx_twoNormOfColumn (
   DenseMtx *mtx,
   int       jcol
) {
double norm;
A2 a2;
A2_setDefaultFields(&a2) ;
DenseMtx_setA2(mtx, &a2);
norm = A2_twoNormOfColumn(&a2, jcol) ;
return (norm); }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   copy column icol of a Dense matrix A to column jcol of a 
   Dense matrix B.  Ai->Bj.

   created -- 98nov12, jwu
   ------------------------------------------------------
*/
void
DenseMtx_colCopy (
   DenseMtx *mtxB,   int jcol, 
   DenseMtx *mtxA,   int icol 
)
{
int nrowA, ncolA, rowincA;
int nrowB, ncolB, rowincB;
double *Ai, *Bj;
int ierr, k;

if (mtxA == NULL || mtxB == NULL) {
  fprintf(stderr,"mtxA and/or mtxB are NULL\n");
  return;
}


if (DENSEMTX_IS_REAL(mtxA) != DENSEMTX_IS_REAL(mtxB)) {
  fprintf(stderr,"mtxA and mtxB do not have the same data type\n");
  return;
}

DenseMtx_dimensions(mtxA, &nrowA, &ncolA);
DenseMtx_dimensions(mtxB, &nrowB, &ncolB);
if (nrowA != nrowB) {
  fprintf(stderr,"Error in DenseMtx_colCopy: nrowA != nrowB\n");
  return;
}
if (icol > ncolA-1 || icol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colCopy: icol is not in [0,ncolA-1]\n");
  return;
}
else if (jcol > ncolB-1 || jcol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colCopy: jcol is not in [0,ncolB-1]\n");
  return;
}

ierr=DenseMtx_column(mtxA, icol, &Ai);
ierr=DenseMtx_column(mtxB, jcol, &Bj);
if (ierr != 1) {
  fprintf(stderr,"Error in DenseMtx_colCopy: ierr!=1 after DenseMtx_column\n");
  return;
}

rowincA=DenseMtx_rowIncrement(mtxA);
rowincB=DenseMtx_rowIncrement(mtxB);

/* copying... */
if (DENSEMTX_IS_COMPLEX(mtxA)) {
   rowincA *= 2;
   rowincB *= 2;
}
for (k=0; k<nrowA; k++) 
  Bj[k*rowincB] = Ai[k*rowincA];

if (DENSEMTX_IS_COMPLEX(mtxA)) {
  Ai+=1;
  Bj+=1;  
  for (k=0; k<nrowA; k++) 
    Bj[k*rowincB] = Ai[k*rowincA];
}

return; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   compute dot product of column icol of a Dense matrix A 
   and column jcol of a Dense matrix B
   prod=Ai^H * Bj.

   created -- 98nov12, jwu
   ------------------------------------------------------
*/

void
DenseMtx_colDotProduct (
   DenseMtx *mtxA,   int icol, 
   DenseMtx *mtxB,   int jcol,
   double *prod
)
{
int nrowA, ncolA, rowincA;
int nrowB, ncolB, rowincB;
double *Ai, *Bj;
int ierr, k;

prod[0]=0.0;
if (!DENSEMTX_IS_REAL(mtxB)) prod[1]=0.0;

if (mtxA == NULL || mtxB == NULL) {
  fprintf(stderr,"mtxA and/or mtxB are NULL\n");
  return;
}

if (DENSEMTX_IS_REAL(mtxA) != DENSEMTX_IS_REAL(mtxB)) {
  fprintf(stderr,"mtxA and mtxB do not have the same data type\n");
  return;
}

DenseMtx_dimensions(mtxA, &nrowA, &ncolA);
DenseMtx_dimensions(mtxB, &nrowB, &ncolB);
if (nrowA != nrowB) {
  fprintf(stderr,"Error in DenseMtx_colDotProduct: nrowA != nrowB\n");
  return;
}
if (icol > ncolA-1 || icol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colDotProduct: "
	  "icol is not in [0,ncolA-1]\n");
  return;
}
else if (jcol > ncolB-1 || jcol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colDotProduct: "
	  "jcol is not in [0,ncolB-1]\n");

  return;
}

ierr=DenseMtx_column(mtxA, icol, &Ai);
ierr=DenseMtx_column(mtxB, jcol, &Bj);
if (ierr != 1) {
  fprintf(stderr,"Error in DenseMtx_colDotProduct: "
	  "ierr!=1 after DenseMtx_column\n");
  return;
}

rowincA=DenseMtx_rowIncrement(mtxA);
rowincB=DenseMtx_rowIncrement(mtxB);

/* inner product */
if (DENSEMTX_IS_COMPLEX(mtxA)) {
   rowincA *= 2;
   rowincB *= 2;
}

for (k=0; k<nrowA; k++) 
  prod[0] += Bj[k*rowincB]*Ai[k*rowincA];

if (DENSEMTX_IS_COMPLEX(mtxA)) {
  Bj+=1;  
  for (k=0; k<nrowA; k++) 
    prod[1] += Bj[k*rowincB]*Ai[k*rowincA];

  Ai+=1;
  for (k=0; k<nrowA; k++) 
    prod[0] += Bj[k*rowincB]*Ai[k*rowincA];

  Bj-=1;  
  for (k=0; k<nrowA; k++) 
    prod[1] -= Bj[k*rowincB]*Ai[k*rowincA];
}

return; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   compute a general axpy with column icol of a Dense matrix A 
   and column jcol of a Dense matrix B.  Ai=alpha*Ai+beta*Bj.

   Note: currently, alpha and beta have to be double, not complex.
 
   created -- 98nov12, jwu
   modified -- 98nov24, jwu
     enable complex alpha and beta scalars
            
   ------------------------------------------------------
*/
void
DenseMtx_colGenAxpy (
   double  *alpha,  DenseMtx *mtxA,   int icol,   
   double  *beta,   DenseMtx *mtxB,   int jcol
)
{
int nrowA, ncolA, rowincA;
int nrowB, ncolB, rowincB;
double *Ai, *Bj, areal, breal, aimag, bimag, tmp;
int ierr, k, j;

if (mtxA == NULL || mtxB == NULL) {
  fprintf(stderr,"mtxA and/or mtxB are NULL\n");
  return;
}

if (DENSEMTX_IS_REAL(mtxA) != DENSEMTX_IS_REAL(mtxB)) {
  fprintf(stderr,"mtxA and mtxB do not have the same data type\n");
  return;
}

DenseMtx_dimensions(mtxA, &nrowA, &ncolA);
DenseMtx_dimensions(mtxB, &nrowB, &ncolB);
if (nrowA != nrowB) {
  fprintf(stderr,"Error in DenseMtx_colGenAxpy: nrowA != nrowB\n");
  return;
}
if (icol > ncolA-1 || icol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colGenAxpy: "
	  "icol is not in [0,ncolA-1]\n");
  return;
}
else if (jcol > ncolB-1 || jcol < 0 ) {
  fprintf(stderr,"Error in DenseMtx_colGenAxpy: "
	  "jcol is not in [0,ncolB-1]\n");
  return;
}

ierr=DenseMtx_column(mtxA, icol, &Ai);
ierr=DenseMtx_column(mtxB, jcol, &Bj);
if (ierr != 1) {
  fprintf(stderr,"Error in DenseMtx_colGenAxpy: "
	  "error after DenseMtx_column\n");
  return;
}

rowincA=DenseMtx_rowIncrement(mtxA);
rowincB=DenseMtx_rowIncrement(mtxB);
areal=*alpha;
breal=*beta;

/* general saxp */
if (DENSEMTX_IS_REAL(mtxA)) {
  if (areal == 0.0)
    for (k=0; k<nrowA; k++) Ai[k*rowincA]=0.0;
  else if (areal != 1.0)
    for (k=0; k<nrowA; k++) Ai[k*rowincA]*=areal;

  if (breal != 0.0) {
    if (breal == 1.0)
      for (k=0; k<nrowA; k++) Ai[k*rowincA]+=Bj[k*rowincB];
    else
      for (k=0; k<nrowA; k++) Ai[k*rowincA]+=breal*Bj[k*rowincB];
  }
}
else {
  rowincA *= 2;
  rowincB *= 2;
  aimag=*(alpha+1);
  bimag=*(beta+1);

  if (areal == 0.0) {
    if (aimag == 0.0) {
      for (k=0; k<nrowA*rowincA; k+=rowincA) {Ai[k]=0.0; Ai[k+1]=0.0;}
    }
    else if (aimag == 1.0){
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; Ai[k]=-Ai[k+1]; Ai[k+1]=tmp;}
    }      
    else {
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; Ai[k]=-aimag*Ai[k+1]; Ai[k+1]=aimag*tmp;}
    }
  }
  else if (areal == 1.0) {
    if (aimag == 1.0) {
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; Ai[k]-=Ai[k+1]; Ai[k+1]+=tmp;}
    }
    else if (aimag != 0.0) {
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; Ai[k]-=aimag*Ai[k+1]; Ai[k+1]+=aimag*tmp;}
    }
  }
  else {
    if (aimag == 0) {
      for (k=0; k<nrowA*rowincA; k+=rowincA) {Ai[k]*=areal; Ai[k+1]*=areal;}
    }
    else if (aimag == 1.0) {
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; Ai[k]=areal*Ai[k]-Ai[k+1]; Ai[k+1]=areal*Ai[k+1]+tmp;}
    }
    else {
      for (k=0; k<nrowA*rowincA; k+=rowincA) 
	{tmp=Ai[k]; 
	 Ai[k]=areal*Ai[k]-aimag*Ai[k+1]; 
	 Ai[k+1]=areal*Ai[k+1]+aimag*tmp;}
    }
  }      

  if (breal == 0.0) {
    if (bimag == 1.0){
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]-=Bj[k+1]; Ai[j+1]+=Bj[k];}
    }      
    else if (bimag != 0.0) {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]-=bimag*Bj[k+1]; Ai[j+1]+=bimag*Bj[k];}
    }
  }
  else if (breal == 1.0) {
    if (bimag == 0.0) {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=Bj[k]; Ai[j+1]+=Bj[k+1];}
    }
    else if (bimag == 1.0) {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=(Bj[k]-Bj[k+1]); Ai[j+1]+=(Bj[k]+Bj[k+1]);}
    }
    else {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=(Bj[k]-bimag*Bj[k+1]); Ai[j+1]+=(Bj[k+1]+bimag*Bj[k]);}
    }
  }
  else {
    if (bimag == 0) {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=Bj[k]*breal; Ai[j+1]+=Bj[k+1]*breal;}
    }
    else if (bimag == 1.0) {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=(breal*Bj[k]-Bj[k+1]); Ai[j+1]+=(breal*Bj[k+1]+Bj[k]);}
    }
    else {
      for (k=j=0; k<nrowB*rowincB; k+=rowincB, j+=rowincA) 
	{Ai[j]+=(breal*Bj[k]-bimag*Bj[k+1]); 
	 Ai[j+1]+=(breal*Bj[k+1]+bimag*Bj[k]);}
    }
  }      
}

return; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   copy col icolA from mtxA into col icolB in mtxB
 
   created -- 98dec10, jwu
   -----------------------------------------------
*/
void
DenseMtx_copyColumn (
   DenseMtx   *mtxB,
   int        icolB,
   DenseMtx   *mtxA,
   int        icolA
) {
double   *colA, *colB ;
int      ii, inc1A, inc1B, iA, iB, nrow ;
/*
   ---------------
   check the input
   ---------------
*/
if (  mtxB == NULL || icolB < 0 || icolB >= mtxB->nrow 
   || mtxA == NULL || icolA < 0 || icolA >= mtxA->nrow 
   || (nrow = mtxA->nrow) != mtxB->nrow ) {
   fprintf(stderr, "\n fatal error in DenseMtx_copyRow(%p,%d,%p,%d)"
           "\n bad input\n", mtxB, icolB, mtxA, icolA) ;
   exit(-1) ;
}
inc1A = mtxA->inc1 ;
inc1B = mtxB->inc1 ;
/*
mtxB->colind[icolB] = mtxA->colind[icolA] ; 
*/
if ( DENSEMTX_IS_REAL(mtxB) && DENSEMTX_IS_REAL(mtxA) ) {
   colA  = mtxA->entries + icolA*mtxA->inc2 ;
   colB  = mtxB->entries + icolB*mtxB->inc2 ;
   for ( ii = iA = iB = 0 ; ii < nrow ; ii++, iA += inc1A, iB += inc1B){
      colB[iB] = colA[iA] ;
   }
} else if ( DENSEMTX_IS_COMPLEX(mtxB) && DENSEMTX_IS_COMPLEX(mtxA) ) {
   colA  = mtxA->entries + 2*icolA*mtxA->inc2 ;
   colB  = mtxB->entries + 2*icolB*mtxB->inc2 ;
   for ( ii = iA = iB = 0 ; ii < nrow ; ii++, iA += inc1A, iB += inc1B){
      colB[2*iB]   = colA[2*iA] ;
      colB[2*iB+1] = colA[2*iA+1] ;
   }
} else {
   fprintf(stderr, "\n fatal error in DenseMtx_copyColumn(%p,%d,%p,%d)"
           "\n mixing real and complex datatype\n", mtxB, icolB, mtxA, icolA);
   exit(-1) ;
}
return ; }
 
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   FrontMtx_solve with column icol of a DenseMtx rhsmtx as 
   right-hand-side.
 
   created -- 98nov24, jwu
   -----------------------------------------------
*/

void
FrontMtx_solveOneColumn (
   FrontMtx        *frontmtx,
   DenseMtx        *solmtx,
   int             jcol,
   DenseMtx        *rhsmtx,
   int             icol,
   SubMtxManager   *mtxmanager,
   double          cpus[],
   int             msglvl,
   FILE            *msgFile
) {

DenseMtx    rhsV, solV;
int         nrowr, nrows, ncol, ierr;

DenseMtx_dimensions (rhsmtx, &nrowr, &ncol);
DenseMtx_dimensions (solmtx, &nrows, &ncol);

/*rhsV = DenseMtx_new() ;*/
DenseMtx_setDefaultFields(&rhsV);
ierr=DenseMtx_initAsSubmatrix (&rhsV, rhsmtx, 0, nrowr-1, icol, icol);
if (ierr < 0) {
   fprintf(stderr, "\n fatal error in FrontMtx_solveOneColumn: ");
   fprintf(stderr, "cannot assign rhsmtx(:,icol) as a sub-DenseMtx");
   exit(-1) ;
}

/*solV = DenseMtx_new() ;*/
DenseMtx_setDefaultFields(&solV);
ierr=DenseMtx_initAsSubmatrix (&solV, solmtx, 0, nrows-1, jcol, jcol);
if (ierr < 0) {
   fprintf(stderr, "\n fatal error in FrontMtx_solveOneColumn: ");
   fprintf(stderr, "cannot assign solmtx(:,jcol) as a sub-DenseMtx");
   exit(-1) ;
}

FrontMtx_solve (frontmtx,&solV,&rhsV,mtxmanager,cpus,msglvl,msgFile);



return; }


/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------
   return the absolute value of a complex number

   created -- 98dec07, ycp
   ----------------------------------------------
*/
 
double 
zabs(double *x)
{
   double v;
   if( fabs(x[0]) > fabs(x[1]) ) {  
      v = fabs(x[1])/fabs(x[0]);
      return ( fabs(x[0])*sqrt(1. + v*v) );
   } 
   else {
     if(x[1] == 0.) 
       return( fabs(x[0]) );
     else {
       v = fabs(x[0])/fabs(x[1]);
       return ( fabs(x[1])*sqrt(1. + v*v) );
     }
   }
}
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   compute the sum of two complex numbers
 
   created -- 98dec07, ycp
   -------------------------------------------
*/

void 
zadd
(double *x, double *y, double *u
) {
   u[0] = x[0] + y[0];
   u[1] = x[1] + y[1];

   return ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   divide two complex numbers                     
 
   created -- 98dec07, ycp 
   ------------------------------------------- 
*/
 
void
zdiv
(double *x, double *y, double *u
) {

   int       flip;
   double    hold, s, t;
   double    x_re, x_im, y_re, y_im;

   x_re = x[0];
   x_im = x[1];
   y_re = y[0];
   y_im = y[1];

   if ( y_im == 0. ){
     u[0] = x_re/y_re;
     u[1] = x_im/y_re;
     return;
   }

   if ( y_re == 0. ){
     u[0] =  x_im/y_im;
     u[1] = -x_re/y_im;
     return;
   }

   flip = 0;
   if ( fabs(y_re) >= fabs(y_im) )
   {
      hold = y_re;
      y_re = y_im;
      y_im = hold;

      hold = x_re;
      x_re = x_im;
      x_im = hold;

      flip = 1;
   }

   s = 1.0/y_re;
   t = 1.0/(y_re+ y_im*(y_im*s) );

   if ( fabs(y_im) >= fabs(s) )
   {
      hold = y_im;
      y_im = s;
      s = hold;
   }

   if ( fabs(x_im) >= fabs(s) )
      u[0] = t*(x_re+ s*(x_im*y_im) );
   else
   {   
      if ( fabs (x_im) >= fabs(y_im) )
         u[0] = t*( x_re + x_im*(s*y_im) );
      else
         u[0] = t*( x_re + y_im*(s*x_im) );
   }
   if ( fabs(x_re) >= fabs(s) )
      u[1] = t*( x_im - s*(x_re*y_im) );
   else
   {
      if ( fabs (x_re) >= fabs(y_im) )
         u[1] = t*(x_im- x_re*(s*y_im) );
      else
         u[1] = t*(x_im- y_im*(s*x_re) );
   }
 
   if ( flip )
      u[1] = -u[1] ;
 
   return ; }
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   multiply two complex numbers
 
   created -- 98dec07, ycp
   -------------------------------------------
*/
 
void
zmul(double *x, double *y,  double *u)
{
double Re, Im;
   Re = x[0]*y[0] - x[1]*y[1] ;
   Im = x[0]*y[1] + x[1]*y[0] ;
   u[0] = Re;
   u[1] = Im;

   return ;
}

 
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   compute the difference of two complex numbers 
 
   created -- 98dec07, ycp
   ---------------------------------------------
*/

void 
zsub
(double *x, double *y, double *u
) {
   u[0] = x[0] - y[0];
   u[1] = x[1] - y[1];

   return ; }

