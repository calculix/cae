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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static char *lakon1,*matname1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*iout1,*nmethod1,
    *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*icmd1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_fn1,calcul_qa1,calcul_cauchy1,*nener1,ikin1,
    *nal=NULL,num_cpus,mt1,*nk1,*ne01,*mortar1,*ielprop1,*kscale1,*neapar1,
    *nebpar1,mscalmethod1=0;

static double *co1,*v1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
    *prestr1,*eme1,*fn1=NULL,*qa1=NULL,*vold1,*veold1,*dtime1,*time1,
    *ttime1,*plicon1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*stiini1,
    *vini1,*ener1,*eei1,*enerini1,*springarea1,*reltime1,
    *thicke1,*emeini1,*prop1,*pslavsurf1,*pmastsurf1,*clearini1,*smscale1,
    enerscal1=0.;

void resultsstr(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *v,double *stn,ITG *inum,double *stx,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *eme,double *emn,
       double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,double *ener,
       double *enern,double *emeini,double *xstaten,double *eei,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime, ITG *ne0,
       double *xforc, ITG *nforc, double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,ITG *nener,
       char *orname,ITG *network,ITG *neapar,ITG *nebpar){

    char *tieset=NULL;
  
    ITG intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
      mt=mi[1]+1,i,j,*itiefac=NULL;

    /*

     calculating the stress integration point values

     iout=2: v is assumed to be known and is used to 
             calculate strains, stresses..., requested results output */
      
    /* variables for multithreading procedure */
    
    ITG sys_cpus,*ithread=NULL;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_SENS");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
	
    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }

// next line is to be inserted in a similar way for all other paralell parts

    if(*ne<num_cpus) num_cpus=*ne;
    
    pthread_t tid[num_cpus];

    /* setting the output variables */

    calcul_fn=0;
    calcul_f=0;
    calcul_qa=0;
    if(iperturb[1]==1){calcul_cauchy=1;}else{calcul_cauchy=0;}

    qa[0]=0.e0;
    qa[1]=0.e0;
    intpointvarm=1;
    ikin=0;

    /* calculating the stresses and material tangent at the 
       integration points; calculating the internal forces */

    if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){

	NNEW(fn1,double,num_cpus*mt**nk);
	NNEW(qa1,double,num_cpus*4);
	NNEW(nal,ITG,num_cpus);

	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;v1=v;
        stx1=stx;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
        nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
        ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
        ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;prestr1=prestr;
        iprestr1=iprestr;eme1=eme;iperturb1=iperturb;iout1=iout;
        vold1=vold;nmethod1=nmethod;veold1=veold;dtime1=dtime;
        time1=time;ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
        plkcon1=plkcon;nplkcon1=nplkcon;xstateini1=xstateini;
        xstiff1=xstiff;xstate1=xstate;npmat1_=npmat_;matname1=matname;
        mi1=mi;ielas1=ielas;icmd1=icmd;ncmat1_=ncmat_;nstate1_=nstate_;
        stiini1=stiini;vini1=vini;ener1=ener;eei1=eei;enerini1=enerini;
        istep1=istep;iinc1=iinc;springarea1=springarea;reltime1=reltime;
        calcul_fn1=calcul_fn;calcul_qa1=calcul_qa;calcul_cauchy1=calcul_cauchy;
        nener1=nener;ikin1=ikin;mt1=mt;nk1=nk;ne01=ne0;thicke1=thicke;
        emeini1=emeini;pslavsurf1=pslavsurf;clearini1=clearini;
        pmastsurf1=pmastsurf;mortar1=mortar;ielprop1=ielprop;prop1=prop;
        kscale1=kscale;neapar1=neapar;
	nebpar1=nebpar;

	/* calculating the stresses */
	
	if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
		printf(" Using up to %" ITGFORMAT " cpu(s) for the stress calculation.\n\n", num_cpus);
	}
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsmechmtstr, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	
	for(i=0;i<mt**nk;i++){
	    fn[i]=fn1[i];
	}
	for(i=0;i<mt**nk;i++){
	    for(j=1;j<num_cpus;j++){
		fn[i]+=fn1[i+j*mt**nk];
	    }
	}
	SFREE(fn1);SFREE(ithread);
	
        /* determine the internal force */

	qa[0]=qa1[0];
	for(j=1;j<num_cpus;j++){
	    qa[0]+=qa1[j*4];
	}

        /* determine the decrease of the time increment in case
           the material routine diverged */

        for(j=0;j<num_cpus;j++){
	    if(qa1[2+j*4]>0.){
		if(qa[2]<0.){
		    qa[2]=qa1[2+j*4];
		}else{
		    if(qa1[2+j*4]<qa[2]){qa[2]=qa1[2+j*4];}
		}
	    }
	}

	SFREE(qa1);
	
	for(j=1;j<num_cpus;j++){
	    nal[0]+=nal[j];
	}

	if(calcul_qa==1){
	    if(nal[0]>0){
		qa[0]/=nal[0];
	    }
	}
	SFREE(nal);
    }

    /* storing results in the .dat file
       extrapolation of integration point values to the nodes
       interpolation of 3d results for 1d/2d elements */

    FORTRAN(resultsprint,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
       stx,ielorien,norien,orab,t1,ithermal,filab,een,iperturb,fn,
       nactdof,iout,vold,nodeboun,ndirboun,nboun,nmethod,ttime,xstate,
       epn,mi,
       nstate_,ener,enern,xstaten,eei,set,nset,istartset,iendset,
       ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
       nelemload,nload,&ikin,ielmat,thicke,eme,emn,rhcon,nrhcon,shcon,
       nshcon,cocon,ncocon,ntmat_,sideload,icfd,inomat,pslavsurf,islavact,
       cdn,mortar,islavnode,nslavnode,ntie,islavsurf,time,ielprop,prop,
       veold,ne0,nmpc,ipompc,nodempc,labmpc,energyini,energy,orname,
       xload,itiefac,pmastsurf,springarea,tieset));
  
  return;

}

/* subroutine for multithreading of resultsmech */

void *resultsmechmtstr(ITG *i){

    ITG indexfn,indexqa,indexnal,nea,neb,nedelta,list,*ilist=NULL;

    indexfn=*i*mt1**nk1;
    indexqa=*i*4;
    indexnal=*i;

    nedelta=(ITG)floor(*ne1/(double)num_cpus);
    
    nea=neapar1[*i]+1;
    neb=nebpar1[*i]+1;
    
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;

    list=0;
    FORTRAN(resultsmech,(co1,kon1,ipkon1,lakon1,ne1,v1,
          stx1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
          ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
          iprestr1,eme1,iperturb1,&fn1[indexfn],iout1,&qa1[indexqa],vold1,
          nmethod1,
          veold1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,nplkcon1,
          xstateini1,xstiff1,xstate1,npmat1_,matname1,mi1,ielas1,icmd1,
          ncmat1_,nstate1_,stiini1,vini1,ener1,eei1,enerini1,istep1,iinc1,
          springarea1,reltime1,&calcul_fn1,&calcul_qa1,&calcul_cauchy1,nener1,
	  &ikin1,&nal[indexnal],ne01,thicke1,emeini1,
	  pslavsurf1,pmastsurf1,mortar1,clearini1,&nea,&neb,ielprop1,prop1,
	  kscale1,&list,ilist,smscale1,&mscalmethod1,&enerscal1));

    return NULL;
}
