/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */

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

#ifdef ARPACK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif
#ifdef MATRIXSTORAGE
   #include "matrixstorage.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif


void randomfieldmain(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,ITG *nmpc,
         ITG *nactdof,ITG *mi,ITG *nodedesi,ITG *ndesi,ITG *istartdesi,
	 ITG *ialdesi,double *co,double *physcon,ITG *isolver,ITG *ntrans,
	 ITG *nk,ITG *inotr,double *trab,char *jobnamec,ITG *nboun,
         double *cs,	 
	 ITG *mcs,ITG *inum,ITG *nmethod,ITG *kode,char *filab,ITG *nstate_,
	 ITG *istep,char *description,char *set,ITG *nset,ITG *iendset,
	 char *output,ITG *istartset,ITG *ialset,double *extnor){
	               
  /* generating the random field for robust assessment */

    char bmat[2]="I", which[3]="SA", howmny[2]="A",*objectset=NULL,*orname=NULL;
  
  ITG nzss,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
      *nx=NULL,*ny=NULL,*nz=NULL,symmetryflag=0,inputformat=0,mei[4],ido,ldz,
      iparam[11],mxiter,info,ncv,jrow,lworkl,nev,ipntr[14],m,
      *select=NULL,rvec=1,i,k,kflag,idesvar,node,j,inorm=0,irand=1,
      iinc=1,mode,noddiam=-1,ngraph,iobject,icoordinate,nrhs=1,ithermal=0;
  
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*au=NULL,*ad=NULL,
         *adb=NULL,*aub=NULL,sigma,*resid=NULL,*workd=NULL,*workl=NULL,tol,
         *temp_array=NULL,*randval=NULL,fmin,fmax,pi,trace,dsum,
	 *d=NULL,*zz=NULL,*acvector=NULL,*randfield=NULL,abserr,
      relerr,*stn=NULL,time=0.;
	         
  //physcon[10]  # of eigenvectors
  //physcon[11]  standarddeviation
  //physcon[12]  correlation length

  /* Assigning values to eigenvalue parameters mei */
  
  mei[1]=*ndesi;  //number of requested Lanczos vectors
  mei[2]=1000;    //number of max iterations
  mei[3]=0;       //matrix storage
  
  /* copying the frequency parameters */

  pi=4.*atan(1.);

  nev=(ITG)(physcon[10]+0.5);
  if(nev>=*ndesi){
      printf("*WARNING in randomfieldmain: number of requested\n");
      printf("         eigenvectors must be smaller than the \n");
      printf("         number of design variables; \n");
      printf("         the number of eigenvectors is set to \n");
      printf("         the number of design variables minus one. \n");
      nev=*ndesi-1;
  }
  ncv=*ndesi;
  //ncv=2*nev;
  mxiter=1000;
  tol=0.01;
  fmin=2*pi*-1;
  fmax=2*pi*-1;
  
  /* prepare for near3d_se; routine "prefilter" can w/o any changes be taken 
     from the routine "filtermain.c" */

  NNEW(xo,double,*ndesi);
  NNEW(yo,double,*ndesi);
  NNEW(zo,double,*ndesi);
  NNEW(x,double,*ndesi);
  NNEW(y,double,*ndesi);
  NNEW(z,double,*ndesi);
  NNEW(nx,ITG,*ndesi);
  NNEW(ny,ITG,*ndesi);
  NNEW(nz,ITG,*ndesi);
  
  for(m=0;m<*ndesi;m++){
      xo[m]=co[3*nodedesi[m]];
      x[m]=xo[m];
      nx[m]=m+1;
      yo[m]=co[3*nodedesi[m]+1];
      y[m]=yo[m];
      ny[m]=m+1;
      zo[m]=co[3*nodedesi[m]+2];
      z[m]=zo[m];
      nz[m]=m+1;
  }
  kflag=2;
  FORTRAN(dsort,(x,nx,ndesi,&kflag));
  FORTRAN(dsort,(y,ny,ndesi,&kflag));
  FORTRAN(dsort,(z,nz,ndesi,&kflag));

  /* determining the structure of the covariance matrix */

  nzss=20000000;
  NNEW(mast1,ITG,nzss);
  NNEW(irows,ITG,1);
  NNEW(icols,ITG,*ndesi);
  NNEW(jqs,ITG,*ndesi+1);
  NNEW(ipointer,ITG,*ndesi);

  mastructrand(icols,jqs,&mast1,&irows,ipointer,&nzss,ndesi,physcon,
               xo,yo,zo,x,y,z,nx,ny,nz);

  SFREE(mast1);SFREE(ipointer);    
  RENEW(irows,ITG,nzss);

  SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
  SFREE(ny);SFREE(nz);
	      
  /* determining the entries of the autocovariance matrix */
  
  NNEW(ad,double,*ndesi);
  NNEW(au,double,nzss);    

  FORTRAN(autocovmatrix,(co,ad,au,jqs,irows,ndesi,nodedesi,physcon));

  /* assign a starting value sigma for the calculation of the eigenvalues
     sigma is set to trace(autocovmatrix) = maximum eigenvalue */
  
  sigma=physcon[11]*physcon[11]**ndesi;
  
  /* Calculate the eigenmodes and eigenvalues of the autocovariance function */

  NNEW(adb,double,*ndesi);
  NNEW(aub,double,nzss); 
  
  for(i=0;i<*ndesi;i++){
     adb[i]=1.0;
  }
    
  /* LU decomposition of the left hand matrix */

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icols,irows,ndesi,&nzss,
                   &symmetryflag,&inputformat,&nzss);
#else
    printf("*ERROR in arpack: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&sigma,icols,irows,ndesi,&nzss,token);
#else
    printf("*ERROR in arpack: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icols,&irows,ndesi,&nzss);
#else
    printf("*ERROR in arpack: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==6){
#ifdef MATRIXSTORAGE
    matrixstorage(ad,&au,adb,aub,&sigma,icols,&irows,ndesi,&nzss,
		  ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
                  lakon,kon,ne,mei,nboun,nmpc,cs,mcs,&ithermal,nmethod);
#else
    printf("*ERROR in arpack: the MATRIXSTORAGE library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icols,irows,ndesi,&nzss,
		   &symmetryflag,&inputformat,jqs,&nzss);
#else
    printf("*ERROR in arpack: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  
  
/* calculating the eigenvalues and eigenmodes */
  
  printf(" Calculating the eigenvalues and the eigenmodes\n\n");

  ido=0;
  ldz=*ndesi;
  iparam[0]=1;
  iparam[2]=mxiter;
  iparam[3]=1;
  iparam[6]=3;
  info=0;

  NNEW(resid,double,*ndesi);
  NNEW(zz,double,(long long)ncv**ndesi);
  NNEW(workd,double,3**ndesi);

  lworkl=ncv*(8+ncv);
  NNEW(workl,double,lworkl);
  FORTRAN(dsaupd,(&ido,bmat,ndesi,which,&nev,&tol,resid,&ncv,zz,&ldz,iparam,ipntr,workd,
      workl,&lworkl,&info));

  NNEW(temp_array,double,*ndesi);
  
  while((ido==-1)||(ido==1)||(ido==2)){
    if(ido==-1){
        FORTRAN(op,(ndesi,&workd[ipntr[0]-1],temp_array,adb,aub,jqs,irows));
    }
    if((ido==-1)||(ido==1)){

      /* solve the linear equation system  */

      if(ido==-1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(temp_array,ndesi);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(temp_array,token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(temp_array,ndesi);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(temp_array,ndesi,&symmetryflag,&nrhs);
#endif
        }
        for(jrow=0;jrow<*ndesi;jrow++){
          workd[ipntr[1]-1+jrow]=temp_array[jrow];
        }
      }
      else if(ido==1){
        if(*isolver==0){
#ifdef SPOOLES
          spooles_solve(&workd[ipntr[2]-1],ndesi);
#endif
        }
        else if(*isolver==4){
#ifdef SGI
          sgi_solve(&workd[ipntr[2]-1],token);
#endif
        }
        else if(*isolver==5){
#ifdef TAUCS
          tau_solve(&workd[ipntr[2]-1],ndesi);
#endif
        }
        else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(&workd[ipntr[2]-1],ndesi,&symmetryflag,&nrhs);
#endif
        }
        for(jrow=0;jrow<*ndesi;jrow++){
          workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
        }
      }

    }

    if(ido==2){
	FORTRAN(op,(ndesi,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
                    adb,aub,jqs,irows));
    }

    FORTRAN(dsaupd,(&ido,bmat,ndesi,which,&nev,&tol,resid,&ncv,zz,&ldz,
                iparam,ipntr,workd,workl,&lworkl,&info));
  }
    
/*--------------------------------------------------------------------*/
/*
   -----------
   free memory
   -----------
*/
  SFREE(temp_array);
  if(*isolver==0){
#ifdef SPOOLES
    spooles_cleanup();
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    sgi_cleanup(token);
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_cleanup();
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
      pardiso_cleanup(ndesi,&symmetryflag);
#endif
  }

  if(info!=0){
    printf("*ERROR in arpack: info=%" ITGFORMAT "\n",info);
    printf("       # of converged eigenvalues=%" ITGFORMAT "\n\n",iparam[4]);
  }         
    
  NNEW(select,ITG,ncv); 
  NNEW(d,double,nev);
  
  FORTRAN(dseupd,(&rvec,howmny,select,d,zz,&ldz,&sigma,bmat,ndesi,which,
      &nev,&tol,resid,&ncv,zz,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
      
      
  /* storing the eigenvectors of the autocovariance matrix in the frd file 
     acvector = autocovariance vector */
  
  NNEW(acvector,double,3**nk);
  
  for(k=nev;k>=1;k--){
      
      /* generating the autocorrelation vector */
      
      for(idesvar=0;idesvar<*ndesi;idesvar++){
	  node=nodedesi[idesvar]-1;
	  for(i=0;i<3;i++){
	      acvector[3*node+i]=zz[(k-1)**ndesi+idesvar]*extnor[3*node+i];
	  }
      }
      
      frd_sen(co,nk,stn,inum,nmethod,kode,filab,&d[k-1],nstate_,
	      istep,
	      &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	      istartset,iendset,ialset,jobnamec,output,
	      acvector,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	      &icoordinate,&inorm,&irand); 
  
  }
  
  /* determining the random values */
  
  NNEW(randval,double,nev);
  FORTRAN(randomval,(randval,&nev));
  
  /* Calculation of the realisation of the random field 
     and the error-measure */
  
  NNEW(randfield,double,*ndesi);

  /* standarddeviation */

  sigma=physcon[11];

  /* generating a random field */

  for(i=0;i<nev;i++){
      for(j=0;j<*ndesi;j++){
	  randfield[j]+=zz[i**ndesi+j]*randval[i]*sqrt(d[i]);
      }
  }

  /* calculating the error measure */

  trace=sigma*sigma**ndesi;
  dsum=0.;
  for(i=0;i<nev;i++){
      dsum+=d[i];
  }
  abserr=trace-dsum;
  relerr=1.-dsum/trace;

  for(k=0;k<*ndesi;k++){
     zz[k]=randfield[k];
  }
      
  /* generating the autocorrelation vector */
  
  for(idesvar=0;idesvar<*ndesi;idesvar++){
      node=nodedesi[idesvar]-1;
      for(i=0;i<3;i++){
	  acvector[3*node+i]=zz[idesvar]*extnor[3*node+i];
      }
  }
      
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&time,nstate_,
	  istep,
	  &iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	  istartset,iendset,ialset,jobnamec,output,
	  acvector,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	  &icoordinate,&inorm,&irand); 
  
  FORTRAN(writerandomfield,(d,&nev,&abserr,&relerr));
  
  SFREE(select);SFREE(workd);SFREE(workl);SFREE(resid);
  SFREE(irows);SFREE(icols);SFREE(jqs);SFREE(au);SFREE(ad);
  SFREE(adb);SFREE(aub);SFREE(d);SFREE(zz);SFREE(randval);
  SFREE(randfield);SFREE(acvector);
       
  return;
  
} 

#endif
