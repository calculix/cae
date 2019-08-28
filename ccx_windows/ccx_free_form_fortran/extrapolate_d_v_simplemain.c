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

static ITG *nface1,*ielfa1,*icyclic1,*ifatie1,*nef1,*num_cpus1;

static double *xrlfa1,*adv1,*advfa1,*hfa1,*c1,*vel1,*volume1;

void extrapolate_d_v_simplemain(ITG *nface,ITG *ielfa,double *xrlfa,double *adv,
				double *advfa,double *hfa,ITG *icyclic,
				double *c,ITG *ifatie,double *vel,ITG *nef,
				double *volume,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */
    
    ielfa1=ielfa,xrlfa1=xrlfa;adv1=adv;advfa1=advfa;hfa1=hfa;
    icyclic1=icyclic;c1=c;ifatie1=ifatie;vel1=vel;nef1=nef;
    volume1=volume;num_cpus1=num_cpus;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)extrapolate_d_v_simple1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of calcgammav1 */

void *extrapolate_d_v_simple1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(extrapolate_d_v_simple,(ielfa1,xrlfa1,adv1,advfa1,hfa1,
				    icyclic1,c1,ifatie1,vel1,nef1,volume1,
         			    &nfacea,&nfaceb));

    return NULL;
}
