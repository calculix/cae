/*     CalculiX - A 3-dimensional finite element program                 */
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

/*     parallel part in resultsini.c	 */


#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *neapar=NULL,*nebpar=NULL,*mt1,*nactdof1;

static double *b1,*v1,*veold1,*accold1,*scal11,*scal21,*cam01,*cam31;

void iniparll(ITG *mt,ITG *nactdof,double *b,double *v,
	      double *veold,double *accold,double *bet,
	      double *gam,double *dtime,double *cam,
	      ITG *nk,ITG *num_cpus){

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    double scal1,scal2;

    scal1=*bet**dtime**dtime;
    scal2=*gam**dtime;
    
    pthread_t tid[*num_cpus];

    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,*num_cpus);
    NNEW(nebpar,ITG,*num_cpus);
    NNEW(cam01,double,*num_cpus);
    NNEW(cam31,double,*num_cpus);

    /* dividing the element number range into num_cpus equal numbers of 
       active entries.  */

    idelta=(ITG)floor(*nk/(double)(*num_cpus));
    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=isum;
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=*nk;
	}
	nebpar[i]=isum;
    }

    /* create threads and wait */
    
    mt1=mt;nactdof1=nactdof;b1=b;v1=v;veold1=veold;accold1=accold;
    scal11=&scal1;scal21=&scal2;
    
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)iniparllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    /* determining the maximum cam[0] */
    
    cam[0]=cam01[0];
    cam[3]=cam31[0];
    for(i=1;i<*num_cpus;i++){
	if(cam01[i]>cam[0]){
	    cam[0]=cam01[i];
	    cam[3]=cam31[i];
	}
    }

    SFREE(ithread);SFREE(neapar);SFREE(nebpar);SFREE(cam01);SFREE(cam31);

}

/* subroutine for multithreading of iniparll */

void *iniparllmt(ITG *i){

    ITG nea,neb,k,j;

    double bnac;

    nea=neapar[*i];
    neb=nebpar[*i];
    
    for(k=nea;k<neb;++k){
	for(j=1;j<*mt1;j++){
	    if(nactdof1[*mt1*k+j]>0){
		bnac=b1[nactdof1[*mt1*k+j]-1];
	    }else{
		continue;
	    }
	    v1[*mt1*k+j]+=*scal11*bnac;
	    if(fabs(*scal11*bnac)>cam01[*i]){
		cam01[*i]=fabs(*scal11*bnac);
		cam31[*i]=nactdof1[*mt1*k+j]-0.5;
	    }
	    veold1[*mt1*k+j]+=*scal21*bnac;
	    accold1[*mt1*k+j]+=bnac;
	}
    }

    return NULL;
}
