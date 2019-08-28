/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

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

ITG *irowpardisoas=NULL,*pointersas=NULL,iparmas[64];
long long ptas[64];
double *aupardisoas=NULL;
ITG nthread_mkl_as=0;

/* factorization, solving and cleaning with PARDISO for
   real unsymmetric matrices */

void pardiso_factor_as(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
                ITG *neq, ITG *nzs, ITG *jq){

  char *env;
  ITG i,j,k,maxfct=1,mnum=1,mtype=11,phase=12,nrhs=1,*perm=NULL,
      msglvl=0,error=0,ifortran,lfortran,index,id;
  long long ndim;
  ITG nthread,nthread_v;
  double *b=NULL,*x=NULL;

  printf(" Factoring the system of radiation equations using the unsymmetric pardiso solver\n");

  iparmas[0]=0;
/* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
   must be done once  */
  if (nthread_mkl_as == 0) {
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
    nthread_mkl_as=nthread;
  }
    
  printf(" number of threads =% d\n\n",nthread_mkl_as);

  for(i=0;i<64;i++){ptas[i]=0;}

  ndim=*neq+(long long)2**nzs;

  NNEW(pointersas,ITG,*neq+1);
  NNEW(irowpardisoas,ITG,ndim);
  NNEW(aupardisoas,double,ndim);

  k=0;

  /* PARDISO requires the matrix to be stored row by row
     aupardisoas contains the entries
     irowpardisoas the corresponding column numbers
         (for each row in ascending order)
     pointersas(i) points to the first entry for row i in 
           field aupardisoas */

  if(*sigma==0.){
    for(i=0;i<*neq;i++){
      pointersas[i]=k+1;
      
      /* lower left triangular matrix */

      for(j=0;j<i;j++){
	ifortran=i+1;
	lfortran=jq[j+1]-jq[j];
	FORTRAN(nident,(&irow[jq[j]-1],&ifortran,&lfortran,&id));
	if(id>0){
	  index=jq[j]+id-2;
	  if(irow[index]==ifortran){
	    irowpardisoas[k]=j+1;
	    aupardisoas[k]=au[index];
	    k++;
	  }
	}
      }

      /* diagonal entry */

      irowpardisoas[k]=i+1;
      aupardisoas[k]=ad[i];
      k++;

      /* upper right triangular matrix */

      for(j=jq[i];j<jq[i+1];j++){
	irowpardisoas[k]=irow[j-1];
	aupardisoas[k]=au[j+*nzs-1];
	k++;
      }
    }
    pointersas[*neq]=k+1;
  }else{
    for(i=0;i<*neq;i++){
      pointersas[i]=k+1;
      
      /* lower left triangular matrix */

      for(j=0;j<i;j++){
	ifortran=i+1;
	lfortran=jq[j+1]-jq[j];
	FORTRAN(nident,(&irow[jq[j]-1],&ifortran,&lfortran,&id));
	if(id>0){
	  index=jq[j]+id-2;
	  if(irow[index]==ifortran){
	    irowpardisoas[k]=j+1;
	    aupardisoas[k]=au[index]-*sigma*aub[index];
	    k++;
	  }
	}
      }

      /* diagonal entry */

      irowpardisoas[k]=i+1;
      aupardisoas[k]=ad[i]-*sigma*adb[i];
      k++;

      /* upper right triangular matrix */

      for(j=jq[i];j<jq[i+1];j++){
	irowpardisoas[k]=irow[j-1];
	aupardisoas[k]=au[j+*nzs-1]-*sigma*aub[j+*nzs-1];
	k++;
      }
    }
    pointersas[*neq]=k+1;
  }

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  return;
}

void pardiso_solve_as(double *b, ITG *neq){

  char *env;
  ITG maxfct=1,mnum=1,mtype=11,phase=33,*perm=NULL,nrhs=1,
      msglvl=0,i,error=0;
  double *x=NULL;

  printf(" Solving the system of radiatio equations using the unsymmetric pardiso solver\n");

  iparmas[0]=0;
/* pardiso_factor has been called befor, MKL_NUM_THREADS=nthread_mkl_as is set*/

  printf(" number of threads =% d\n\n",nthread_mkl_as);

  NNEW(x,double,*neq);

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  for(i=0;i<*neq;i++){b[i]=x[i];}
  SFREE(x);

  return;
}

void pardiso_cleanup_as(ITG *neq){

  ITG maxfct=1,mnum=1,mtype=11,phase=-1,*perm=NULL,nrhs=1,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  SFREE(irowpardisoas);
  SFREE(aupardisoas);
  SFREE(pointersas);

  return;
}

void pardiso_main_as(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, ITG *icol, ITG *irow, 
         ITG *neq, ITG *nzs, ITG *jq){

  if(*neq==0) return;

  pardiso_factor_as(ad,au,adb,aub,sigma,icol,irow, 
             neq,nzs,jq);

  pardiso_solve_as(b,neq);

  pardiso_cleanup_as(neq);

  return;
}

#endif

