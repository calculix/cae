/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                     */

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
#include <pthread.h>
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
#ifdef PASTIX
#include "pastix.h"
#endif


void projectgradmain(ITG *nobject,char **objectsetp,double **dgdxglobp,
		     double *g0,ITG *ndesi,ITG *nodedesi,ITG *nk,ITG *isolver,
		     double *co, double *xdesi, double *distmin, ITG *nconstraint){
               
  /* generating the projected gradient vector for constraint 
     optimization */

  char *objectset=NULL;
  
  ITG nzss,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
    symmetryflag=0,inputformat=0,i,iconst,iter,nactiveold,*ipoactiold=NULL,
    *iconstactiold=NULL,*ipoacti=NULL,*iconstacti=NULL,nactive=0,nnlconst,
    iscaleflag,nrhs=1,*inameacti=NULL,*inameactiold=NULL,iobject,node;
             
  double *au=NULL,*ad=NULL,*adb=NULL,*aub=NULL,sigma=0,*rhs=NULL,
    *vector=NULL,*xlambd=NULL,*xtf=NULL,*objnorm=NULL,*dgdxglob=NULL,
    *dgdxglobcpy=NULL;   
  
  objectset=*objectsetp;dgdxglob=*dgdxglobp;
  
  /* check if any constraints are active */

  /* in the field "ipoacti" and the variable "nnlconst" the 
     characteristics of the N-matrix are stored:
     -nnlconst is the number of nonlinear constraints
     -in the first nnlconst entries of ipoacti the position of the 
     nonlinear constraints (e.g. mass) is inserted
     -in the following entries the position of the linear 
     constraints (wall thickness constraints and fixation of nodes) 
     in nodedesi(i) is stored (i=ipoacti(j))
     -the field "iconstacti" contains the information whether the  
     constraint is LE or GE 
     > -1 is LE
     >  1 is GE
     -in the field "inameacti" the name of the active constraint is 
     saved via a pointer to the field objectset */
  
  NNEW(objnorm,double,*nobject);
  NNEW(ipoacti,ITG,*nobject**ndesi);
  NNEW(inameacti,ITG,*nobject**ndesi);  
  NNEW(iconstacti,ITG,*nobject**ndesi);

  /* estimate nactive on the basis of the function values of the 
     constraints */
     
  if(nconstraint>0){
  
    FORTRAN(checkconstraint,(nobject,objectset,g0,&nactive,&nnlconst,
			     ipoacti,ndesi,dgdxglob,nk,nodedesi,iconstacti,objnorm,
			     inameacti));	
  }

  RENEW(ipoacti,ITG,nactive);
  RENEW(inameacti,ITG,nactive);
  RENEW(iconstacti,ITG,nactive);
     
  if(nactive>0){
     
    iscaleflag=1;     
    FORTRAN(scalesen,(dgdxglob,nobject,nk,nodedesi,ndesi,
		      objectset,&iscaleflag,ipoacti,&nnlconst,&nactive));
      
    *nobject=*nobject+1; 
    RENEW(dgdxglob,double,2**nk**nobject);
    RENEW(objectset,char,324**nobject);
    for(i=324*(*nobject-1);i<324**nobject;i++)objectset[i]=' ';
     
    nactiveold=nactive+1;
    iter=0;
    nzss=20000000;
    NNEW(mast1,ITG,nzss);
    NNEW(irows,ITG,1);
    NNEW(icols,ITG,nactive);
    NNEW(jqs,ITG,nactive+1);
    NNEW(ad,double,nactive);
    NNEW(au,double,nzss);    
    NNEW(adb,double,nactive);
    NNEW(aub,double,nzss); 
    NNEW(rhs,double,nactive);
    NNEW(xlambd,double,nactive);
    NNEW(xtf,double,nactive);
    NNEW(vector,double,*ndesi);
    NNEW(ipoactiold,ITG,nactive);
    NNEW(iconstactiold,ITG,nactive);

    while((nactive<nactiveold)&&(nactive>0)){
  
      nactiveold=nactive;
      iter=iter+1;
     
      /* determining the structure of the N-Matrix */

      nzss=20000000;
      RENEW(mast1,ITG,nzss);
      RENEW(irows,ITG,1);
      RENEW(icols,ITG,nactive);
      RENEW(jqs,ITG,nactive+1);
      NNEW(ipointer,ITG,nactive);
     
      mastructnmatrix(icols,jqs,&mast1,&irows,ipointer,&nzss,&nactive,&nnlconst);

      RENEW(irows,ITG,nzss);
      SFREE(ipointer);
     
      /* determining the entries of the N-Matrix */
  
      RENEW(ad,double,nactive);
      RENEW(au,double,nzss);    
  
      FORTRAN(nmatrix,(ad,au,jqs,irows,ndesi,nodedesi,dgdxglob,&nactive,
		       nobject,&nnlconst,ipoacti,nk));

      /* Calculate inverse of the N-matrix */

      RENEW(adb,double,nactive);
      RENEW(aub,double,nzss); 

      for(i=0;i<nactive;i++){
	adb[i]=1.0;
      }
       
      /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
	spooles_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
		       &symmetryflag,&inputformat,&nzss);
#else
	printf("*ERROR in projectgrad: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,token);
#else
	printf("*ERROR in projectgrad: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(ad,&au,adb,aub,&sigma,icols,&irows,&nactive,&nzss);
#else
	printf("*ERROR in projectgrad: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
		       &symmetryflag,&inputformat,jqs,&nzss);
#else
	printf("*ERROR in projectgrad: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(ad,au,adb,aub,&sigma,icols,irows,&nactive,&nzss,
			   &symmetryflag,&inputformat,jqs,&nzss);
#else
	printf("*ERROR in projectgrad: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
	
      /* solve the system nactive-times */

      for(i=0;i<nactive;i++){
	xtf[i]=0.00;
      }
	
      RENEW(rhs,double,nactive);
      RENEW(xlambd,double,nactive);
      RENEW(xtf,double,nactive);
      RENEW(vector,double,*ndesi);
	
      FORTRAN(preprojectgrad,(vector,ndesi,nodedesi,dgdxglob,&nactive,
			      nobject,&nnlconst,ipoacti,nk,rhs,&iconst,objectset,xtf));
	
      /* Calculate the projected gradient and the lagrange multipliers */
	
      for(iconst=1;iconst<=nactive;iconst++){
	    
	for(i=0;i<nactive;i++){
	  rhs[i]=0.00;
	}
	    
	rhs[iconst-1]=1.0;
	    
	/* solve the system */
	    
	if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(rhs,&nactive);
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	  sgi_solve(rhs,token);
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(rhs,&nactive);
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  pardiso_solve(rhs,&nactive,&symmetryflag,&inputformat,&nrhs);
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  pastix_solve(rhs,&nactive,&symmetryflag,&nrhs);
#endif
	}
	    
	for(i=0;i<*ndesi;i++){
	  vector[i]=0.00;
	}
	    
	/* carry out matrix multiplications */
	    
	FORTRAN(projectgrad,(vector,ndesi,nodedesi,dgdxglob,&nactive,
			     nobject,&nnlconst,ipoacti,nk,rhs,&iconst,objectset,xlambd,xtf,
			     objnorm));
	    
      }
	
      /* clean the system */
	
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
	pardiso_cleanup(&nactive,&symmetryflag,&inputformat);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
#endif
      }
	
      /* write the results of the iteration in the dat-file */
	
      FORTRAN(writelm,(&iter,xlambd,&nactive,&nnlconst,objectset,nobject,
		       ipoacti,iconstacti,inameacti));
	
      /* check if the langrange multipliers of the active constraints
	 have the "correct" sign 
	 ("correct" in the sense of active constraints) */
	
      RENEW(ipoactiold,ITG,nactive);
      RENEW(iconstactiold,ITG,nactive);
      RENEW(inameactiold,ITG,nactive);
	
      FORTRAN(checkprojectgrad,(&nactiveold,&nactive,ipoacti,ipoactiold,
				objectset,xlambd,&nnlconst,iconstacti,
				iconstactiold,inameacti,inameactiold,g0,
				nobject));
	
    }

    SFREE(inameactiold);
    
    /* normalization of projected gradient and documentation of results */
     
    if(nactive==0){
	 
      *nobject=*nobject-1;  
      RENEW(dgdxglob,double,2**nk**nobject);
      RENEW(objectset,char,324**nobject);
     
    }
   
  }
  
  FORTRAN(postprojectgrad,(ndesi,nodedesi,dgdxglob,&nactive,nobject,
			   &nnlconst,ipoacti,nk,&iconst,objectset,iconstacti,inameacti));
                
  /* Filtering of sensitivities a second time to account for dependency
     between the control field and the designvariables*/

  NNEW(dgdxglobcpy,double,2**nk**nobject);
  memcpy(&dgdxglobcpy[0],&dgdxglob[0],sizeof(double)*2**nk**nobject);
  
  iscaleflag=2;	     	     
  FORTRAN(scalesen,(dgdxglobcpy,nobject,nk,nodedesi,ndesi,
		    objectset,&iscaleflag,ipoacti,&nnlconst,&nactive));
     
  filtermain(co,dgdxglobcpy,nobject,nk,nodedesi,ndesi,objectset,
             xdesi,distmin);

  for(iobject=0;iobject<*nobject;iobject++){
    for(i=0;i<*ndesi;i++){
      node=nodedesi[i]-1;
      dgdxglob[2*node+1+2*iobject**nk]=dgdxglobcpy[2*node+1+2*iobject**nk];
    }
  }
       
  iscaleflag=3;
  FORTRAN(scalesen,(dgdxglob,nobject,nk,nodedesi,ndesi,
		    objectset,&iscaleflag,ipoacti,&nnlconst,&nactive));
  
  SFREE(mast1);SFREE(irows);SFREE(icols);SFREE(jqs);
  SFREE(ad);SFREE(au);SFREE(adb);SFREE(aub);SFREE(rhs);SFREE(xlambd);
  SFREE(xtf);SFREE(vector);SFREE(ipoactiold);SFREE(iconstactiold);
  SFREE(objnorm);SFREE(ipoacti);SFREE(iconstacti);SFREE(inameacti);
  SFREE(dgdxglobcpy);
  
  *objectsetp=objectset;*dgdxglobp=dgdxglob;
  
  return;
  
} 

#endif
