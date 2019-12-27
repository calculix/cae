/*  DenseMtx_mmm.c */

#include "../Iter.h"
 
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   performs the matrix-matrix operations
   C =  beta*C + alpha*(A)*(B)
   A, B and C must be column major.

   Parameters ---
      A_opt -- form of op( A ) to be used
               = 'N' or 'n',  op( A ) = A.
               = 'T' or 't',  op( A ) = A'.
               = 'C' or 'c',  op( A ) = A*.
      B_opt -- form of op( B ) to be used
               = 'N' or 'n',  op( B ) = B.
               = 'T' or 't',  op( B ) = B'.
               = 'C' or 'c',  op( B ) = B*.
   return values ---
      1 -- normal return
     -1 -- C, A, B, alpha or beta is NULL
     -2 -- type of A, B or C are invalid
     -3 -- row or column of A and B is not match
     -4 -- invalid option for A or B

   created -- 98dec11, ycp
   -------------------------------------------
*/

int DenseMtx_mmm(
  char     *A_opt,
  char     *B_opt,
  double   *beta,
  DenseMtx *mtxC, 
  double   *alpha,
  DenseMtx *mtxA, 
  DenseMtx *mtxB
  ) 
{

int nrowA, ncolA, rowincA, colincA;
int nrowB, ncolB, rowincB, colincB;
int nrowC, ncolC, rowincC, colincC;
int ierr, i, k, j, l;
double *Ai, *Bj, *Ci, r_alpha, r_beta, r_temp, im_temp, im_alpha, im_beta;
double  one[2]={1.0, 0.0}, zero[2]={0.0, 0.0}, aconj[2], bconj[2] ;
double temp[2]={0.0, 0.0}, result[2]={1.0, 0.0} ;

if ( beta == NULL || alpha == NULL || mtxC == NULL ||
     mtxA == NULL || mtxB  == NULL ){ 
  fprintf(stderr, "\n fatal error in Input"
          "\n one or more of beta, alpha, mtxC, mtxB and"
          " mtxA is NULL\n") ;
  return(-1) ;
}
if ( (DENSEMTX_IS_REAL(mtxA) != DENSEMTX_IS_REAL(mtxB)) ||
     (DENSEMTX_IS_REAL(mtxA) != DENSEMTX_IS_REAL(mtxC)) ){
  fprintf(stderr,"mtxA, mtxB and mtxC do not have the same data type\n");
  return(-2);
}
 
DenseMtx_dimensions(mtxA, &nrowA, &ncolA);
DenseMtx_dimensions(mtxB, &nrowB, &ncolB);
DenseMtx_dimensions(mtxC, &nrowC, &ncolC);

rowincA=DenseMtx_rowIncrement(mtxA);
colincA=DenseMtx_columnIncrement(mtxA);
rowincB=DenseMtx_rowIncrement(mtxB);
colincB=DenseMtx_columnIncrement(mtxB);
rowincC=DenseMtx_rowIncrement(mtxC);
colincC=DenseMtx_columnIncrement(mtxC);

r_alpha=*alpha;
r_beta =*beta;
r_temp =*temp;

if ( B_opt[0] == 'N' || B_opt[0] == 'n' ){        
  if (A_opt[0] == 'N' || A_opt[0] == 'n'){/*Form C := beta*c+alpha*A*B*/
    if (ncolA != nrowB || nrowC != nrowA || ncolC != ncolB) {
      fprintf(stderr,"Error in Input DenseMtx_mmm\n");
      return(-3);
    }
  } else if ( (A_opt[0] == 'T' || A_opt[0] == 't') ||
              (A_opt[0] == 'C' || A_opt[0] == 'c')  ){
    if (nrowA != nrowB || nrowC != ncolA || ncolC != ncolB) {
       fprintf(stderr,"Error in Input DenseMtx_mmm\n");
       exit(-3);
    }
  } else {
    fprintf(stderr,"Invalid option for mtxA\n");
    return(-4);
  }
} else if ( (B_opt[0] == 'T' || B_opt[0] == 't') ||
            (B_opt[0] == 'C' || B_opt[0] == 'c')  ){
  if (A_opt[0] == 'N' || A_opt[0] == 'n'){
    if (ncolA != ncolB || nrowC != nrowA || ncolC != nrowB) {
      fprintf(stderr,"Error in Input DenseMtx_mmm\n");
      return(-3);
    }
  } else if ( (A_opt[0] == 'T' || A_opt[0] == 't') ||
              (A_opt[0] == 'C' || A_opt[0] == 'c')  ){
    if (nrowA != ncolB || nrowC != ncolA || ncolC != nrowB) {
      fprintf(stderr,"Error in Input DenseMtx_mmm\n");
      return(-3);
    }
  } else {
    fprintf(stderr,"Invalid option for mtxA\n");
    return(-4);
  }
} else {
  fprintf(stderr,"Invalid option for mtxB\n");
  return(-4);
}

if (DENSEMTX_IS_REAL(mtxA)) {
  if ( r_alpha == *zero ) {
    if( r_beta == *zero ) {
       DenseMtx_zero (mtxC);
    } else {
       DenseMtx_scale(mtxC,&r_beta);
    }
    return(1);
  }
  if ( B_opt[0] == 'N' || B_opt[0] == 'n' ){        
    if (A_opt[0] == 'N' || A_opt[0] == 'n'){/*Form C := beta*c+alpha*A*B*/
      for (i=0; i<nrowA; i++){
        ierr=DenseMtx_row(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<ncolB; j++){
          ierr=DenseMtx_column(mtxB, j, &Bj);
          r_temp = 0.0;
          for (k=0; k<ncolA; k++){
            r_temp += Ai[k*colincA]*Bj[k*rowincB];
          }
          if( r_beta == *zero ){
            Ci[j*colincC] = r_alpha*r_temp; 
          } else {
            Ci[j*colincC] = r_alpha*r_temp + r_beta*Ci[j*colincC]; 
          }
        }
      }
    } else {/* Form  C := alpha*AT*B + beta*C. */
      for (i=0; i<ncolA; i++){
        ierr=DenseMtx_column(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<ncolB; j++){
          ierr=DenseMtx_column(mtxB, j, &Bj);
          r_temp = 0.0;
          for (k=0; k<nrowA; k++){
            r_temp += Ai[k*rowincA]*Bj[k*rowincB];
          }
          if( r_beta == *zero ){
            Ci[j*colincC] = r_alpha*r_temp;
          } else {
            Ci[j*colincC] = r_alpha*r_temp + r_beta*Ci[j*colincC];
          }
        }
      }
    }
  } else {
    if (A_opt[0] == 'N' || A_opt[0] == 'n'){/*Form  C := alpha*A*B'+beta*C */
      for (i=0; i<nrowA; i++){
        ierr=DenseMtx_row(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          r_temp = 0.0;
          for (k=0; k<ncolA; k++){
            r_temp += Ai[k*colincA]*Bj[k*colincB];
          }
          if( r_beta == *zero ){
            Ci[j*colincC] = r_alpha*r_temp;
          } else {
            Ci[j*colincC] = r_alpha*r_temp + r_beta*Ci[j*colincC];
          }
        }
      }
    } else { /* Form  C := alpha*A'*B' + beta*C */
      for (i=0; i<ncolA; i++){
        ierr=DenseMtx_column(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          r_temp = 0.0;
          for (k=0; k<nrowA; k++){
            r_temp += Ai[k*rowincA]*Bj[k*colincB];
          }
          if( r_beta == *zero ){
            Ci[j*colincC] = r_alpha*r_temp;
          } else {
            Ci[j*colincC] = r_alpha*r_temp + r_beta*Ci[j*colincC];
          }
        }
      }
    }
  }
} else { /* complex case */
  rowincA *= 2;
  rowincB *= 2;
  rowincC *= 2;
  colincA *= 2;
  colincB *= 2;
  colincC *= 2;
  im_alpha=*(alpha+1);
  im_beta =*(beta+1);
  im_temp =*(temp+1);

  if ( r_alpha == *zero && im_alpha == *zero ) {
    if( r_beta == *zero && im_beta == *zero ) {
       DenseMtx_zero (mtxC);
    } else {
       DenseMtx_scale(mtxC,beta);
    }
    return(1);
  }
  if ( B_opt[0] == 'N' || B_opt[0] == 'n' ){
    if (A_opt[0] == 'N' || A_opt[0] == 'n'){/*Form C := beta*c+alpha*A*B*/
      for (i=0; i<nrowA; i++){
        ierr=DenseMtx_row(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<ncolB; j++){
          temp[0] = 0.0;
          temp[1] = 0.0;
          ierr=DenseMtx_column(mtxB, j, &Bj);
          for (k=0; k<ncolA; k++){
            zmul(&Ai[k*colincA],&Bj[k*rowincB],result);
            zadd(temp,result,temp);
          }
          if( r_beta == *zero && im_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zmul(alpha,temp,temp);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    } else {/* Form  C := alpha*AT*B + beta*C. */
      for (i=0; i<ncolA; i++){
        ierr=DenseMtx_column(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<ncolB; j++){
          ierr=DenseMtx_column(mtxB, j, &Bj);
          temp[0] = 0.0;
          temp[1] = 0.0;
          for (k=0; k<nrowA; k++){
            if (A_opt[0] == 'C' || A_opt[0] == 'c'){
              /* Form  C := alpha*conjg( A')*B + beta*C. */
              aconj[0] = Ai[k*rowincA];
              aconj[1] = -Ai[k*rowincA+1];
              zmul(aconj,&Bj[k*rowincB],result);
            } else {
              zmul(&Ai[k*rowincA],&Bj[k*rowincB],result);
            }
            zadd(temp,result,temp);
          }
          if( r_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(alpha,temp,temp);
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    }
  } else if ( B_opt[0] == 'T' || B_opt[0] == 'T' ){
    if (A_opt[0] == 'N' || A_opt[0] == 'n'){/*Form  C := alpha*A*B'+beta*C */
      for (i=0; i<nrowA; i++){
        ierr=DenseMtx_row(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          temp[0] = 0.0;
          temp[1] = 0.0;
          for (k=0; k<ncolA; k++){
            zmul(&Ai[k*colincA],&Bj[k*colincB],result);
            zadd(temp,result,temp);
          }
          if( r_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(alpha,temp,temp);
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    } else { /* Form  C := alpha*A'*B' + beta*C */
      for (i=0; i<ncolA; i++){
        ierr=DenseMtx_column(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          temp[0] = 0.0;
          temp[1] = 0.0;
          for (k=0; k<nrowA; k++){
            if (A_opt[0] == 'C' || A_opt[0] == 'c'){
               /* Form  C := alpha*conjg( A')*B' + beta*C. */

              aconj[0] = Ai[k*rowincA];
              aconj[1] = (-1)*Ai[k*rowincA+1];
              zmul(aconj,&Bj[k*rowincB],result);
            } else {
              zmul(&Ai[k*rowincA],&Bj[k*rowincB],result);
            }
            zadd(temp,result,temp);
          }
          if( r_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(alpha,temp,temp);
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    }
  } else {
    if (A_opt[0] == 'N' || A_opt[0] == 'n'){
                      /*Form  C := alpha*A*conjg(B')+beta*C */
      for (i=0; i<nrowA; i++){
        ierr=DenseMtx_row(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          temp[0] = 0.0;
          temp[1] = 0.0;
          for (k=0; k<ncolA; k++){
            bconj[0] = Bj[k*colincB];
            bconj[1] = -Bj[k*colincB+1];
            zmul(&Ai[k*colincA],bconj,result);
            zadd(temp,result,temp);
          }
          if( r_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(alpha,temp,temp);
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    } else { /* Form  C := alpha*A'*conjg(B') + beta*C */
      for (i=0; i<ncolA; i++){
        ierr=DenseMtx_column(mtxA, i, &Ai);
        ierr=DenseMtx_row(mtxC, i, &Ci);
        for (j=0; j<nrowB; j++){
          ierr=DenseMtx_row(mtxB, j, &Bj);
          temp[0] = 0.0;
          temp[1] = 0.0;
          for (k=0; k<nrowA; k++){
            bconj[0] =  Bj[k*colincB];
            bconj[1] = -Bj[k*colincB+1];
            if (A_opt[0] == 'C' || A_opt[0] == 'c'){
               /* Form  C := alpha*conjg( A')*conjg(B') + beta*C. */
              aconj[0] = Ai[k*rowincA];
              aconj[1] = -Ai[k*rowincA+1];
              zmul(aconj,bconj,result);
            } else {
              zmul(&Ai[k*rowincA],bconj,result);
            }
            zadd(temp,result,temp);
          }
          if( r_beta == *zero ){
            zmul(alpha,temp,&Ci[j*colincC]);
          } else {
            zmul(alpha,temp,temp);
            zmul(beta,&Ci[j*colincC],&Ci[j*colincC]);
            zadd(temp,&Ci[j*colincC],&Ci[j*colincC]);
          }
        }
      }
    }
  }
}
return(1); }
