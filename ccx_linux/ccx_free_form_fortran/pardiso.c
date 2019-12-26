/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2019 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "pardiso.h"

ITG *icolpardiso=NULL,*pointers=NULL,iparm[64];
long long pt[64];
double *aupardiso=NULL;
/* double dparm[64];  not used */
ITG nthread_mkl=0;
/* char envMKL[32];   moved to pardiso.h */

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){

  char *env;
/*  char env1[32]; */
  ITG i,j,k,l,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
      msglvl=0,error=0,*irowpardiso=NULL,kflag,kstart,n,ifortran,
      lfortran,index,id,k2;
  ITG ndim,nthread,nthread_v;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0){
      printf(" Factoring the system of equations using the symmetric pardiso solver\n");
  }else{
      printf(" Factoring the system of equations using the unsymmetric pardiso solver\n");
  }

  iparm[0]=0;
  iparm[1]=3;
/* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
   must be done once  */
  if (nthread_mkl == 0) {
    nthread=1;
    env=getenv("MKL_NUM_THREADS");
    if(env) {
      nthread=atoi(env);}
    else {
      env=getenv("OMP_NUM_THREADS");
      if(env) {nthread=atoi(env);}
    }
    env=getenv("CCX_NPROC_EQUATION_SOLVER");
    if(env) {
      nthread_v=atoi(env);
      if (nthread_v <= nthread) {nthread=nthread_v;}
    }
    if (nthread < 1) {nthread=1;}
    sprintf(envMKL,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);  
    putenv(envMKL);
    nthread_mkl=nthread;
  }
    
  printf(" number of threads =% d\n\n",nthread_mkl);

  for(i=0;i<64;i++){pt[i]=0;}

  if(*symmetryflag==0){

      /* symmetric matrix; the subdiagonal entries are stored
         column by column in au, the diagonal entries in ad;
         pardiso needs the entries row per row */      

      mtype=-2;
      
      ndim=*neq+*nzs;
      
      NNEW(pointers,ITG,*neq+1);
      NNEW(icolpardiso,ITG,ndim);
      NNEW(aupardiso,double,ndim);
      
      k=ndim;
      l=*nzs;
      
      if(*sigma==0.){
	  pointers[*neq]=ndim+1;
	  for(i=*neq-1;i>=0;--i){
	      for(j=0;j<icol[i];++j){
		  icolpardiso[--k]=irow[--l];
		  aupardiso[k]=au[l];
	      }
	      pointers[i]=k--;
	      icolpardiso[k]=i+1;
	      aupardiso[k]=ad[i];
	  }
      }
      else{
	  pointers[*neq]=ndim+1;
	  for(i=*neq-1;i>=0;--i){
	      for(j=0;j<icol[i];++j){
		  icolpardiso[--k]=irow[--l];
		  aupardiso[k]=au[l]-*sigma*aub[l];
	      }
	      pointers[i]=k--;
	      icolpardiso[k]=i+1;
	      aupardiso[k]=ad[i]-*sigma*adb[i];
	  }
      }
  }else{
//      mtype=11;
      mtype=1;

      if(*inputformat==3){

          /* off-diagonal terms  are stored column per
             column from top to bottom in au;
             diagonal terms are stored in ad  */

	  ndim=*neq+*nzs;
	  NNEW(pointers,ITG,*neq+1);
	  NNEW(irowpardiso,ITG,ndim);	  
	  NNEW(icolpardiso,ITG,ndim);
	  NNEW(aupardiso,double,ndim);
	  
	  k=0;
	  k2=0;
	  for(i=0;i<*neq;i++){
	      for(j=0;j<icol[i];j++){
		  if(au[k]>1.e-12||au[k]<-1.e-12){
		      icolpardiso[k2]=i+1;
		      irowpardiso[k2]=irow[k];
		      aupardiso[k2]=au[k];
		      k2++;		  
		  }
		  k++;	      
	      }	  
	  }  
	  /* diagonal terms */  
	  for(i=0;i<*neq;i++){
	      icolpardiso[k2]=i+1;
	      irowpardiso[k2]=i+1;
	      aupardiso[k2]=ad[i];
	      k2++;	  
	  }
	  ndim=k2;
	  
	  /* pardiso needs the entries row per row; so sorting is
	     needed */ 
	  
	  kflag=2;
	  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			    &ndim,&kflag));
	  
	  /* sorting each row */
	  
	  k=0;
	  pointers[0]=1;
	  for(i=0;i<*neq;i++){
	      j=i+1;
	      kstart=k;
	      do{
		  if(irowpardiso[k]!=j ){
		      n=k-kstart;		  
		      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				       &n,&kflag));
		      pointers[i+1]=k+1;
		      break;  
		  }else{
		      if(k+1==ndim){
			  n=k-kstart+1;	  
			  FORTRAN(isortid,(&icolpardiso[kstart],
                                  &aupardiso[kstart],&n,&kflag));
			  break;	       
		      }else{
			  k++;	       
		      }  
		  }
	      }while(1);
	  }
	  pointers[*neq]=ndim+1;
	  SFREE(irowpardiso);

      }else if(*inputformat==1){
	  
          /* lower triangular matrix is stored column by column in
             au, followed by the upper triangular matrix row by row;
             the diagonal terms are stored in ad */

          /* reordering lower triangular matrix */

	  ndim=*nzs;
	  NNEW(pointers,ITG,*neq+1);
	  NNEW(irowpardiso,ITG,ndim);
	  NNEW(icolpardiso,ITG,ndim);
	  NNEW(aupardiso,double,ndim);
	  
	  k=0;
	  for(i=0;i<*neq;i++){
	      for(j=0;j<icol[i];j++){
		  icolpardiso[k]=i+1;
		  irowpardiso[k]=irow[k];
		  aupardiso[k]=au[k];
		  k++;
	      }
	  }
	  
	  /* pardiso needs the entries row per row; so sorting is
	     needed */
	  
	  kflag=2;
	  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
	  &ndim,&kflag));
	  
	  /* sorting each row */
	  
	  k=0;
	  pointers[0]=1;
	  if(ndim>0){
	      for(i=0;i<*neq;i++){
		  j=i+1;
		  kstart=k;
		  do{

	  	      /* end of row reached */

		      if(irowpardiso[k]!=j){
			  n=k-kstart;
			  FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					   &n,&kflag));
			  pointers[i+1]=k+1;
			  break;
		      }else{

		          /* end of last row reached */

			  if(k+1==ndim){
			      n=k-kstart+1;
			      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					       &n,&kflag));
			      break;
			  }else{

			      /* end of row not yet reached */

			      k++;
			  }
		      }
		  }while(1);
	      }
	  }
	  pointers[*neq]=ndim+1;
	  SFREE(irowpardiso);

	  /* composing the matrix: lower triangle + diagonal + upper triangle */

	  ndim=*neq+2**nzs;
	  RENEW(icolpardiso,ITG,ndim);
	  RENEW(aupardiso,double,ndim);
	  k=ndim;
	  for(i=*neq-1;i>=0;i--){
	      l=k+1;
	      for(j=jq[i+1]-1;j>=jq[i];j--){
		  icolpardiso[--k]=irow[j-1];
		  aupardiso[k]=au[j+*nzs3-1];
	      }
	      icolpardiso[--k]=i+1;
	      aupardiso[k]=ad[i];
	      for(j=pointers[i+1]-1;j>=pointers[i];j--){
		  icolpardiso[--k]=icolpardiso[j-1];
		  aupardiso[k]=aupardiso[j-1];
	      }
	      pointers[i+1]=l;
	  }
	  pointers[0]=1;
      }
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  return;
}

void pardiso_solve(double *b, ITG *neq,ITG *symmetryflag,ITG *nrhs){

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
      msglvl=0,i,error=0;
  double *x=NULL;

  if(*symmetryflag==0){
      printf(" Solving the system of equations using the symmetric pardiso solver\n");
  }else{
      printf(" Solving the system of equations using the unsymmetric pardiso solver\n");
  }

  if(*symmetryflag==0){
      mtype=-2;
  }else{
//      mtype=11;
      mtype=1;
  }
//  iparm[0]=0;
  iparm[1]=3;
/* pardiso_factor has been called befor, MKL_NUM_THREADS=nthread_mkl is set*/

  printf(" number of threads =% d\n\n",nthread_mkl);

  NNEW(x,double,*nrhs**neq);

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,nrhs,iparm,&msglvl,
                   b,x,&error));

  for(i=0;i<*nrhs**neq;i++){b[i]=x[i];}
  SFREE(x);

  return;
}

void pardiso_cleanup(ITG *neq,ITG *symmetryflag){

  ITG maxfct=1,mnum=1,phase=-1,*perm=NULL,nrhs=1,mtype,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0){
      mtype=-2;
  }else{
//      mtype=11;
      mtype=1;
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  SFREE(icolpardiso);
  SFREE(aupardiso);
  SFREE(pointers);

  return;
}

void pardiso_main(double *ad, double *au, double *adb, double *aub, 
         double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	 ITG *jq, ITG *nzs3,ITG *nrhs){

  if(*neq==0) return;

  pardiso_factor(ad,au,adb,aub,sigma,icol,irow, 
		 neq,nzs,symmetryflag,inputformat,jq,nzs3);

  pardiso_solve(b,neq,symmetryflag,nrhs);

  pardiso_cleanup(neq,symmetryflag);

  return;
}

#endif

