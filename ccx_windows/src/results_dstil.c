/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                          */

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
#include "mortar.h"

static char *lakon1,*matname1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*iout1,*nmethod1,
    *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*icmd1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_fn1,calcul_qa1,calcul_cauchy1,*nener1,ikin1,
    *nal=NULL,
    num_cpus,mt1,*nk1,*ne01,*mortar1,
    *ielprop1,*kscale1,
    *neapar=NULL,*nebpar=NULL,*irowtloc1,*jqtloc1,*islavelinv1;

static double *co1,*v1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
    *prestr1,*eme1,*fn1=NULL,*qa1=NULL,*vold1,*veold1,*dtime1,*time1,
    *ttime1,*plicon1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*stiini1,
    *vini1,*ener1,*eei1,*enerini1,*springarea1,*reltime1,
    *thicke1,*emeini1,*prop1,
    *pslavsurf1,*pmastsurf1,*clearini1,
    *autloc1;
/**
 *  function to calculate internal forces with transformed basis functions
 * Author: Saskia Sitzmann
 * @see: results.c
 *
 *  [in] pslavsurf        field storing  position xil, etal and weight for integration point on slave side
 *  [in] pmastsurf         field storing position and etal for integration points on master side 
 *  [in] islavact                (i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
 *  [in] cdn                ?
 *  [in] islavnode        field storing the nodes of the slave surface 
 *  [in] nslavnode        (i)pointer into field isalvnode for contact tie i
 *  [in] islavelinv       (i)==0 if there is no slave node in the element, >0 otherwise
 *  [in] autloc                transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$  
 *  [in] irowtloc                field containing row numbers of autloc
 *  [in] jqtloc                pointer into field irowtloc 
 *  [in] nboun2                  number of transformed SPCs
 *  [in] ndirboun2        (i) direction of transformed SPC i 
 *  [in] nodeboun2               (i) node of transformed SPC i
 *  [in] xboun2                  (i) value of transformed SPC i
 *  [in] nmpc2                number of transformed mpcs
 *  [in] ipompc2                 (i) pointer to nodempc and coeffmpc for transformed MPC i
 *  [in] nodempc2                nodes and directions of transformed MPCs
 *  [in] coefmpc2                coefficients of transformed MPCs
 *  [in] labmpc2                 transformed mpc labels
 *  [in] ikboun2                 sorted dofs idof=8*(node-1)+dir for transformed SPCs
 *  [in] ilboun2                 transformed SPC numbers for sorted dofs
 *  [in] ikmpc2                 sorted dofs idof=8*(node-1)+dir for transformed MPCs
 *  [in] ilmpc2                transformed SPC numbers for sorted dofs
 **/

void results_dstil(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *v,double *stn,ITG *inum,double *stx,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,double *t1,ITG *ithermal,double *prestr,ITG *iprestr,
       char *filab,double *eme,double *emn,
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
       double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
       char *orname,ITG *network,ITG *ipobody,double *xbody,ITG *ibody,
       char *typeboun,
       ITG *islavelinv,double *autloc,ITG *irowtloc,ITG *jqtloc,
       ITG *nboun2,ITG *ndirboun2,ITG *nodeboun2,double *xboun2,
       ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,char *labmpc2,
       ITG *ikboun2,ITG *ilboun2,ITG *ikmpc2,ITG *ilmpc2){

    ITG intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
        intpointvart,mt=mi[1]+1,i,j;

    /*

     calculating integration point values (strains, stresses,
     heat fluxes, material tangent matrices and nodal forces)

     storing the nodal and integration point results in the
     .dat file

     iout=-2: v is assumed to be known and is used to
              calculate strains, stresses..., no result output
              corresponds to iout=-1 with in addition the
              calculation of the internal energy density
     iout=-1: v is assumed to be known and is used to
              calculate strains, stresses..., no result output;
              is used to take changes in SPC's and MPC's at the
              start of a new increment or iteration into account
     iout=0: v is calculated from the system solution
             and strains, stresses.. are calculated, no result output
     iout=1:  v is calculated from the system solution and strains,
              stresses.. are calculated, requested results output
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

    envloc = getenv("CCX_NPROC_RESULTS");
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

    if(*ne<num_cpus) num_cpus=*ne;
    
    pthread_t tid[num_cpus];
    
    /* 1. nodewise storage of the primary variables
       2. determination which derived variables have to be calculated */

    resultsini(nk,v,ithermal,filab,iperturb,f,fn,
       nactdof,iout,qa,vold,b,nodeboun2,ndirboun2,
       xboun2,nboun2,ipompc2,nodempc2,coefmpc2,labmpc2,nmpc2,nmethod,cam,neq,
       veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
       &ikin,&intpointvart,typeboun,&num_cpus);

    /* calculating the stresses and material tangent at the 
       integration points; calculating the internal forces */

    if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){

        /* determining the element bounds in each thread */

	NNEW(neapar,ITG,num_cpus);
	NNEW(nebpar,ITG,num_cpus);
	elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

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
	kscale1=kscale;islavelinv1=islavelinv;
        autloc1=autloc;jqtloc1=jqtloc;irowtloc1=irowtloc;

        /* calculating the stresses */
	
	printf(" Using up to %" ITGFORMAT " cpu(s) for the stress calculation.\n\n", num_cpus);
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsmechmt_dstil, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)
	  pthread_join(tid[i], NULL);
	
	for(i=0;i<mt**nk;i++){
	    fn[i]=fn1[i];
	}
	for(i=0;i<mt**nk;i++){
	    for(j=1;j<num_cpus;j++){
		fn[i]+=fn1[i+j*mt**nk];
	    }
	}
	SFREE(fn1);SFREE(ithread);SFREE(neapar);SFREE(nebpar);
	
        /* determine the internal force */

	qa[0]=qa1[0];
	for(j=1;j<num_cpus;j++){
	    qa[0]+=qa1[j*4];
	}

        /* determine the decrease of the time increment in case
           the material routine diverged */

	qa[2]=qa1[2];
        for(j=1;j<num_cpus;j++){
	    if(qa1[2+j*4]>0.){
		if(qa[2]<0.){
		    qa[2]=qa1[2+j*4];
		}else{
		    if(qa1[2+j*4]<qa[2]){
		      qa[2]=qa1[2+j*4];}
		}
	    }
	}

        /* maximum change in creep strain increment in the
           present time increment */

	qa[3]=qa1[3];
        for(j=1;j<num_cpus;j++){
	    if(qa1[3+j*4]>0.){
		if(qa[3]<0.){
		    qa[3]=qa1[3+j*4];
		}else{
		    if(qa1[3+j*4]>qa[3]){
		      qa[3]=qa1[3+j*4];}
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

    /* calculating the matrix system internal force vector */

    resultsforc(nk,f,fn,nactdof,ipompc2,nodempc2,
		coefmpc2,labmpc2,nmpc2,mi,fmpc,&calcul_fn,&calcul_f,
		&num_cpus);

  return;

}

/* subroutine for multithreading of resultsmech */

void *resultsmechmt_dstil(ITG *i){

    ITG indexfn,indexqa,indexnal,nea,neb,list1,*ilist1=NULL;

    indexfn=*i*mt1**nk1;
    indexqa=*i*4;
    indexnal=*i;

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

    list1=0;     
    FORTRAN(resultsmech_dstil,(co1,kon1,ipkon1,lakon1,ne1,v1,
          stx1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
          ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
          iprestr1,eme1,iperturb1,&fn1[indexfn],iout1,&qa1[indexqa],vold1,
          nmethod1,
          veold1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,nplkcon1,
          xstateini1,xstiff1,xstate1,npmat1_,matname1,mi1,ielas1,icmd1,
          ncmat1_,nstate1_,stiini1,vini1,ener1,eei1,enerini1,istep1,iinc1,
          springarea1,reltime1,&calcul_fn1,&calcul_qa1,&calcul_cauchy1,nener1,
	  &ikin1,&nal[indexnal],ne01,thicke1,emeini1,pslavsurf1,pmastsurf1,
	  mortar1,clearini1,&nea,&neb,ielprop1,prop1,
	  kscale1,&list1,ilist1,
          islavelinv1,autloc1,irowtloc1,jqtloc1));

    return NULL;
}
