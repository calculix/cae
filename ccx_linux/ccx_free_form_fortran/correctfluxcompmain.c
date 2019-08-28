/*     CalculiX - A 3-dimensional finite element program                 */
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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *nef1,*ipnei1,*neifa1,*neiel1,*ielfa1,*ifabou1,*ielmatf1,
    *mi1,*ntmat1_,*num_cpus1;

static double *flux1,*vfa1,*advfa1,*area1,*vel1,*xlet1,*xle1,*shcon1;

void correctfluxcompmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,double *flux,
			 double *vfa,double *advfa,double *area,double *vel,
			 double *xlet,ITG *ielfa,double *xle,ITG *ifabou,
			 ITG *ielmatf,ITG *mi,double *shcon,ITG *ntmat_,
			 ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */

    nef1=nef;ipnei1=ipnei;neifa1=neifa;neiel1=neiel;flux1=flux;vfa1=vfa;
    advfa1=advfa;area1=area;vel1=vel;xlet1=xlet;ielfa1=ielfa;xle1=xle;
    ifabou1=ifabou;ielmatf1=ielmatf;mi1=mi;shcon1=shcon;ntmat1_=ntmat_;
    num_cpus1=num_cpus;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)correctfluxcomp1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of materialdatacfd1 */

void *correctfluxcomp1mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(correctfluxcomp,(nef1,ipnei1,neifa1,neiel1,flux1,vfa1,advfa1,
			      area1,vel1,xlet1,ielfa1,xle1,ifabou1,ielmatf1,
			      mi1,shcon1,ntmat1_,&nefa,&nefb));

    return NULL;
}
