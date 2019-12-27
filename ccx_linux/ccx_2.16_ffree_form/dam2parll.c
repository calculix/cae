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

/*     A parallel copy of arrays which depent on active element 	 */


#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *neapar=NULL,*nebpar=NULL;

static double *dampwk1,*cv1,*cvini1,*aux21;

void dam2parll(double *dampwk,double *cv,double *cvini,
		    double *aux2,ITG *neq0,ITG *num_cpus){

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    pthread_t tid[*num_cpus];

    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,*num_cpus);
    NNEW(nebpar,ITG,*num_cpus);
    NNEW(dampwk1,double,*num_cpus);

    /* dividing the element number range into num_cpus equal numbers of 
       active entries.  */

    idelta=(ITG)floor(*neq0/(double)(*num_cpus));
    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=isum;
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=*neq0;
	}
	nebpar[i]=isum;
    }

    /* create threads and wait */
    
    cv1=cv;cvini1=cvini;aux21=aux2;
    
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)dam2parllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    for(i=0;i<*num_cpus;i++){
	(*dampwk)+=dampwk1[i];
    }

    SFREE(ithread);SFREE(neapar);SFREE(nebpar);SFREE(dampwk1);

}

/* subroutine for multithreading of copyparll */

void *dam2parllmt(ITG *i){

    ITG nea,neb,k;

    nea=neapar[*i];
    neb=nebpar[*i];
    
    for(k=nea;k<neb;++k){
	dampwk1[*i]+=-(cv1[k]+cvini1[k])*aux21[k]/2.;
    }

    return NULL;
}
