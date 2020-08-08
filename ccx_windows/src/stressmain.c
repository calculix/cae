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
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"

static char *lakon1,*matname1,*filabl1,*labmpc1,*set1,*prlab1,*prset1,
  *sideload1,*orname1,*objectset1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
  *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*nmethod1,
  *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*ncmat1_,*nstate1_,
  *istep1,*iinc1,*istartdesi1,*ialdesi1,
  num_cpusd,*ne01,*mortar1,*ielprop1,*ndesi1,*nodedesi1,idesvar1,
  *nobject1,iobject1,*nk1,*nactdof1,*nodeboun1,*ndirboun1,
  *nboun1,*ipompc1,*nodempc1,*nmpc1,*neq1,idesvara1,idesvarbref1,
  *ikboun1,*ilboun1,*ncocon1,*nset1,*istartset1,*iendset1,*ialset1,
  *nprint1,*inotr1,*ntrans1,*nelemload1,*nload1,*ikmpc1,*ilmpc1,*nforc1,
  *nshcon1,*icfd1,*inomat1,*islavact1,*islavnode1,*nslavnode1,*ntie1,
  *islavsurf1,kscale1,network1,nestart1,neend1,*jqs1,*irows1,*nasym1,
  *isolver1,nodeset1,num_cpuse,*ialdesi1,
  *ialnk1,*ialeneigh1,*ialnneigh1,*ipos1,*nodedesired1,*istarteneigh1,
  *istartnneigh1,*nactdofred1,*nactdofinv1,*mt1,*istartnk1,
  *iponod2dto3d1,*ipobody1,*ibody1,*nbody1;
    
    
static double *co1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
  *prestr1,*vold1,*veold1,*dtime1,*time1,*xdesi1,
  *ttime1,*plicon1,*plkcon1,*xstateini1,*xstate1,*stiini1,
  *vini1,*enerini1,*springarea1,*reltime1,*thicke1,*emeini1,
  *prop1,*pslavsurf1,*pmastsurf1,*clearini1,*distmin1,*g01,*dgdx1,
  *een1,*f1,*xboun1,
  *coefmpc1,*cam1,*accold1,*bet1,*gam1,*epn1,*enern1,*xstaten1,*cocon1,*qfx1,
  *qfn1,*trab1,*fmpc1,*xforc1,*shcon1,*xload1,*xloadold1,*cdn1,*energyini1,
  *energy1,*emn1,*stn1,*b1,*dgdxdy1,*dv1,*dstn1,*dstx1,*expks1,*dgdu1,
  dispmin1,*conew1,*xbody1;
    

void stressmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
		double *v,double *stn,ITG *inum,double *stx,double *elcon,
		ITG *nelcon,double *rhcon,ITG *nrhcon,double *alcon,
		ITG *nalcon,double *alzero,ITG *ielmat,ITG *ielorien,
		ITG *norien,double *orab,ITG *ntmat_,double *t0,double *t1,
		ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
		double *eme,double *emn,double *een,ITG *iperturb,double *f,
		double *fn,ITG *nactdof,ITG *iout,double *qa,double *vold,
		ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
		ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
		ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
		double *accold,double *bet,double *gam,double *dtime,
		double *time,double *ttime,double *plicon,ITG *nplicon,
		double *plkcon,ITG *nplkcon,double *xstateini,double *xstiff,
		double *xstate,ITG *npmat_,double *epn,char *matname,ITG *mi,
		ITG *ielas,ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
		double *vini,ITG *ikboun,ITG *ilboun,double *ener,
		double *enern,double *emeini,double *xstaten,double *eei,
		double *enerini,double *cocon,ITG *ncocon,char *set,ITG *nset,
		ITG *istartset,ITG *iendset,ITG *ialset,ITG *nprint,
		char *prlab,char *prset,double *qfx,double *qfn,double *trab,
		ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
		ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,double *springarea,
		double *reltime,ITG *ne0,double *xforc,ITG *nforc,
		double *thicke,double *shcon,ITG *nshcon,char *sideload,
		double *xload,double *xloadold,ITG *icfd,ITG *inomat,
		double *pslavsurf,double *pmastsurf,ITG *mortar,ITG *islavact,
		double *cdn,ITG *islavnode,ITG *nslavnode,ITG *ntie,
		double *clearini,ITG *islavsurf,ITG *ielprop,double *prop,
		double *energyini,double *energy,double *distmin,ITG *ndesi,
		ITG *nodedesi,ITG *nobject,char *objectset,double *g0,
		double *dgdx,double *sti,double *df,ITG *nactdofinv,ITG *jqs,
		ITG *irows,ITG *idisplacement,ITG *nzs,char *jobnamec,
		ITG *isolver,ITG *icol,ITG *irow,ITG *jq,ITG *kode,double *cs,
		char *output,ITG *istartdesi,ITG *ialdesi,double *xdesi,
		char *orname,ITG *icoordinate,ITG *iev,double *d,double *z,
		double *au,double *ad,double *aub,double*adb,
		ITG *cyclicsymmetry,ITG *nzss,ITG *nev,ITG *ishapeenergy,
		double *fint,ITG *nlabel,ITG *igreen,ITG *nasym,ITG *iponoel,
		ITG *inoel,ITG *nodedesiinv,double *dgdxglob,double *df2,
		double *dgdxdy,ITG *nkon,ITG *iponod2dto3d,ITG *iponk2dto3d,
		ITG *ics,ITG *mcs,ITG *mpcend,ITG *noddiam,double *duds,
		ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody){

  char *filabl=NULL;

  ITG nener=0,i,j,k,m,iobject,symmetryflag=0,inputformat=0,
    mt=mi[1]+1,idesvar,nea,neb,nodeset,lmax,
    kscale=1,network=0,*neinset=NULL,
    nepar,isum,idelta,*neapar=NULL,*nebpar=NULL,nestart,neend,num_cpus,
    l,idesvarbref,idesvara,idesvarb,inode,node,idof,nrhs=1,
    *istartnk,*ialnk,*istartnneigh=NULL,*ialnneigh=NULL,*ichecknodes=NULL,
    *icheckelems=NULL,*istarteneigh=NULL,*ialeneigh=NULL,neielemtot,
    *nkinsetinv=NULL,ndesired,*nodedesired=NULL,neqred,
    *nactdofred=NULL,ipos;
	
  double sigma=0.,*dstn=NULL,*b=NULL,dispmin=1.e-8,
    expks,*dgdu=NULL,*dv=NULL,*dstx=NULL,*conew=NULL;	
    
  if(*nasym!=0){symmetryflag=2;inputformat=1;}

  /* variables for multithreading procedure */
    
  ITG sys_cpus,*ithread=NULL;
  char *env,*envloc,*envsys;
    
  num_cpus=0;
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
    
  pthread_t tid[num_cpus];

  /* factor the system */

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,&nzs[2]);
#else
    printf("*ERROR in objectivemain_se: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
    printf("*ERROR in objectivemain_se: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
    printf("*ERROR in objectivemain_se: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
    printf("*ERROR in objectivemain_se: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
	    
  /* loop over all objective functions */

  for(m=0;m<*nobject;m++){
    if(strcmp1(&objectset[m*324],"STRESS")==0){
      iobject=m+1;

      NNEW(filabl,char,87**nlabel);
      for(i=0;i<87**nlabel;i++){strcpy1(&filabl[i]," ",1);}
      strcpy1(&filabl[174],"S   ",4);
	    
      /* deactivating all elements which are not part of
	 the target function */
	   
      NNEW(neinset,ITG,*ne);
      NNEW(nkinsetinv,ITG,*nk);

      FORTRAN(actideactistr,(set,nset,istartset,iendset,ialset,objectset,
			     ipkon,&iobject,ne,neinset,iponoel,inoel,&nepar,
			     nkinsetinv,nk));

      /* storing the elements to which each node belongs
	 in field ialnk */

      NNEW(istartnk,ITG,*nk+1);
      NNEW(ialnk,ITG,*nkon);

      FORTRAN(createialnk,(nk,iponoel,inoel,istartnk,ialnk,ipkon));

      RENEW(ialnk,ITG,istartnk[*nk]-1);
	
      /* storing the nodes of the neighboring elements of node nk 
	 in field ialnneigh */
	
      NNEW(istartnneigh,ITG,*nk+1);
      NNEW(ialnneigh,ITG,20*(istartnk[*nk]-1));
      NNEW(ichecknodes,ITG,*nk);
	
      FORTRAN(createnodeneigh,(nk,istartnk,ialnk,
			       istartnneigh,ialnneigh,ichecknodes,lakon,ipkon,kon,
			       nkinsetinv,&neielemtot));	    
	    
      RENEW(ialnneigh,ITG,istartnneigh[*nk]-1); 
      SFREE(ichecknodes);  
      SFREE(nkinsetinv); 

      NNEW(istarteneigh,ITG,*nk+1);
      NNEW(ialeneigh,ITG,neielemtot);
      NNEW(icheckelems,ITG,*ne);

      FORTRAN(createelemneigh,(nk,iponoel,inoel,istartnneigh,
			       ialnneigh,icheckelems,istarteneigh,ialeneigh));
	
      RENEW(ialeneigh,ITG,istarteneigh[*nk]-1);
      SFREE(icheckelems);   
 	                  	    
      /* determining the nodal bounds in each thread */

      if(nepar<num_cpus){num_cpuse=nepar;}else{num_cpuse=num_cpus;}

      NNEW(neapar,ITG,num_cpuse);
      NNEW(nebpar,ITG,num_cpuse);
	    
      idelta=nepar/num_cpuse;
	    
      /* dividing the range from 1 to the number of active elements */
	    
      isum=0;
      for(i=0;i<num_cpuse;i++){
	neapar[i]=neinset[isum]-1;
	if(i!=num_cpuse-1){
	  isum+=idelta;
	}else{
	  isum=nepar;
	}
	nebpar[i]=neinset[isum-1]-1;
      }

      /* FORTRAN convention */

      nestart=neapar[0]+1;
      neend=nebpar[num_cpuse-1]+1;

      SFREE(neinset);
	    
      /* OBJECTIVE: STRESS */

      /* calculating the stress in the unperturbed state */
  
      NNEW(v,double,mt**nk);
      NNEW(fn,double,mt**nk);
      NNEW(stn,double,6**nk);
      NNEW(inum,ITG,*nk);
      NNEW(stx,double,6*mi[0]**ne);
      NNEW(eei,double,6*mi[0]**ne);
	    
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      *iout=2;
      *icmd=3;
	    
      resultsstr(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		 ielorien,norien,orab,ntmat_,t0,t1,ithermal,
		 prestr,iprestr,filabl,eme,emn,een,iperturb,
		 f,fn,nactdof,iout,qa,vold,b,nodeboun,
		 ndirboun,xboun,nboun,ipompc,
		 nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		 bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
		 xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
		 ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		 xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		 ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		 nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
		 reltime,ne0,xforc,nforc,thicke,shcon,nshcon,
		 sideload,xload,xloadold,icfd,inomat,pslavsurf,pmastsurf,
		 mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		 islavsurf,ielprop,prop,energyini,energy,&kscale,
		 &nener,orname,&network,neapar,nebpar,ipobody,ibody,xbody,
		 nbody);
	    
      *icmd=0;
	    
      SFREE(v);SFREE(fn);SFREE(eei);
	    

      /* check for the existence of a target node set */

      /* calculating the objective function */

      nodeset=0;
      for(i=0;i<*nset;i++){
	if(strcmp1(&objectset[m*324+162]," ")==0) continue;
	if(strcmp2(&objectset[m*324+162],&set[i*81],81)==0){
	  nodeset=i+1;
	  break;
	}
      }
      FORTRAN(objective_stress,(&nodeset,istartset,iendset,
				ialset,nk,&idesvar,&iobject,mi,g0,
				nobject,stn,objectset,&expks));

      /* reduce the size of the loops */
	    
      NNEW(nodedesired,ITG,*ndesi);
      ndesired=0;
	    
      for(idesvar=1;idesvar<=*ndesi;idesvar++){
              
	node=nodedesi[idesvar-1];
	nea=istarteneigh[node-1];
	neb=istarteneigh[node]-1;
	
	if(neb>=nea){
	        
	  nodedesired[ndesired]=idesvar;
	  ndesired=ndesired+1;
	              
	}
	
      }
	    
      RENEW(nodedesired,ITG,ndesired);
      NNEW(nactdofred,ITG,*neq);
      neqred=0;
	    
      for(idof=0;idof<*neq;idof++){

	inode=nactdofinv[idof];
	//	idir=inode-mt*(inode/mt);
	node=inode/mt+1;
              
	nea=istarteneigh[node-1];
	neb=istarteneigh[node]-1;
	    
	if(neb>=nea){
	        
	  nactdofred[neqred]=idof;
	  neqred=neqred+1;
	      
	}
	
      }
	    
      RENEW(nactdofred,ITG,neqred);	       

      NNEW(dgdu,double,*neq);
      NNEW(dv,double,mt**nk*num_cpus);
      NNEW(dstn,double,6**nk*num_cpus);
      NNEW(dstx,double,6*mi[0]**ne*num_cpus);
      NNEW(conew,double,3**nk*num_cpus);
	    
      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
      stn1=stn;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
      nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
      ielmat1=ielmat;ielorien1=ielorien;norien1=norien;
      orab1=orab;ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;
      prestr1=prestr;iprestr1=iprestr;filabl1=filabl;
      iperturb1=iperturb;vold1=vold;nmethod1=nmethod;dtime1=dtime;
      time1=time;ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
      plkcon1=plkcon;nplkcon1=nplkcon;xstateini1=xstateini;
      xstate1=xstate;npmat1_=npmat_;matname1=matname;mi1=mi;
      ielas1=ielas;ncmat1_=ncmat_;nstate1_=nstate_;stiini1=stiini;
      vini1=vini;emeini1=emeini;enerini1=enerini;istep1=istep;	  
      iinc1=iinc;springarea1=springarea;reltime1=reltime;ne01=ne0;
      thicke1=thicke;pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;
      mortar1=mortar;clearini1=clearini;ielprop1=ielprop;
      prop1=prop;kscale1=kscale;iobject1=iobject;objectset1=objectset;
      g01=g0;dgdx1=dgdx;nasym1=nasym;distmin1=distmin;idesvar1=idesvar;
      dv1=dv;dstn1=dstn;dstx1=dstx;ialdesi1=ialdesi;ialnk1=ialnk;
      dgdu1=dgdu;ialeneigh1=ialeneigh;ialnneigh1=ialnneigh;expks1=&expks; 	    
      stx1=stx;dispmin1=dispmin;ipos1=&ipos;istartnneigh1=istartnneigh;	    
      istarteneigh1=istarteneigh;conew1=conew;nodedesired1=nodedesired;
      nodedesi1=nodedesi;istartdesi1=istartdesi;xdesi1=xdesi;
      nactdofred1=nactdofred;nactdofinv1=nactdofinv;mt1=&mt;
      istartnk1=istartnk;ndesi1=ndesi;iponod2dto3d1=iponod2dto3d;

      /* Variation of the coordinates of the designvariables */
	    
      printf(" original number of design variables: %" ITGFORMAT " \n", *ndesi);
      printf(" reduced number of design variables: %" ITGFORMAT " \n\n", ndesired);
	    
      printf(" Using up to %" ITGFORMAT " cpu(s) for the sensitivity of the stress w.r.t. the coordinates.\n\n", num_cpus);

      lmax=ndesired/num_cpus;

      /* deviding the design variables in sets of
	 num_cpus variables */

      for(l=0;l<lmax+1;l++){
	if(l<lmax){
	  num_cpusd=num_cpus;
	}else{
	  num_cpusd=ndesired-lmax*num_cpus;
	  if(num_cpusd==0){break;}
	}

	ipos=l*num_cpus;
	
	NNEW(ithread,ITG,num_cpusd);
	
	for(i=0;i<num_cpusd;i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)stress_sen_dx_robustmt, (void *)&ithread[i]);
	}
	
	for(i=0;i<num_cpusd;i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);

      }
	    
      /* Variation of the displacements */

      printf(" original number of dofs: %" ITGFORMAT " \n", *neq);
      printf(" reduced number of dofs: %" ITGFORMAT " \n\n", neqred);
	    
      printf(" Using up to %" ITGFORMAT " cpu(s) for the sensitivity of the stress w.r.t. the displacements.\n\n", num_cpus);
	    	    
      lmax=neqred/num_cpus;

      /* deviding the design variables in sets of
	 num_cpus variables */

      for(l=0;l<lmax+1;l++){
	if(l<lmax){
	  num_cpusd=num_cpus;
	}else{
	  num_cpusd=neqred-lmax*num_cpus;
	  if(num_cpusd==0){break;}
	}

	ipos=l*num_cpus;
	
	NNEW(ithread,ITG,num_cpusd);
	
	for(i=0;i<num_cpusd;i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)stress_sen_dv_robustmt, (void *)&ithread[i]);
	}
	
	for(i=0;i<num_cpusd;i++)  pthread_join(tid[i], NULL);
	
	SFREE(ithread);

      }
	    
      SFREE(b);SFREE(dv);SFREE(dstn);SFREE(dstx);SFREE(stx);
      SFREE(conew);SFREE(istartnneigh);SFREE(ialnneigh);SFREE(ialnk);
      SFREE(istartnk);SFREE(ialeneigh);SFREE(istarteneigh);
      SFREE(nactdofred);SFREE(nodedesired);

      //	    /* Multiplication of dg/du with K^-1 */		    
      //	    
      //	    /* solve the system */
      //	
      //	      if(*isolver==0){
      //#ifdef SPOOLES
      //	          spooles_solve(dgdu,&neq[1]);
      //#endif
      //	      }
      //	      else if(*isolver==4){
      //#ifdef SGI
      //	          sgi_solve(dgdu,token);
      //#endif
      //	      }
      //	      else if(*isolver==5){
      //#ifdef TAUCS
      //	          tau_solve(dgdu,&neq[1]);
      //#endif
      //	      }
      //	      else if(*isolver==7){
      //#ifdef PARDISO
      //	          pardiso_solve(dgdu,&neq[1],&symmetryflag,&nrhs);
      //#endif
      //	      }
	      
      /* calculation of total differential */
	    
      /*      FORTRAN(objective_stress_tot_robust,(dgdx,df,ndesi,&iobject,jqs,
	      irows,dgdu,duds,neq));*/
            					    
      /* reactivating all elements */

      for(i=0;i<*ne;i++){
	if(ipkon[i]<-1) ipkon[i]=-2-ipkon[i];
      }

      SFREE(inum);SFREE(stn);SFREE(filabl);
      SFREE(neapar);SFREE(nebpar);SFREE(dgdu);
	    
    }else if(strcmp1(&objectset[m*324],"STR2ND")==0){
      iobject=m+1;

      NNEW(filabl,char,87**nlabel);
      for(i=0;i<87**nlabel;i++){strcpy1(&filabl[i]," ",1);}
      strcpy1(&filabl[174],"S   ",4);
	    
      /* deactivating all elements which are not part of
	 the target function */
	   
      NNEW(neinset,ITG,*ne);

      FORTRAN(actideactistr,(set,nset,istartset,iendset,ialset,objectset,
			     ipkon,&iobject,ne,neinset,iponoel,inoel,&nepar,
			     nkinsetinv,nk));

      /* determining the nodal bounds in each thread */

      if(nepar<num_cpus){num_cpuse=nepar;}else{num_cpuse=num_cpus;}

      NNEW(neapar,ITG,num_cpuse);
      NNEW(nebpar,ITG,num_cpuse);
	    
      idelta=nepar/num_cpuse;
	    
      /* dividing the range from 1 to the number of active elements */
	    
      isum=0;
      for(i=0;i<num_cpuse;i++){
	neapar[i]=neinset[isum]-1;
	if(i!=num_cpuse-1){
	  isum+=idelta;
	}else{
	  isum=nepar;
	}
	nebpar[i]=neinset[isum-1]-1;
      }
	    
      /* translating the bounds of the ranges to real node numbers */

      for(i=0;i<num_cpuse;i++){
	printf("i=%d,neapar=%d,nebpar=%d\n",i,neapar[i],nebpar[i]);
      }

      /* FORTRAN convention */

      nestart=neapar[0]+1;
      neend=nebpar[num_cpuse-1]+1;

      SFREE(neinset);
	    
      /* OBJECTIVE: STRESS 2ND ORDER TREATMENT */

      /* calculating the stress in the unperturbed state */
  
      NNEW(v,double,mt**nk);
      NNEW(fn,double,mt**nk);
      NNEW(stn,double,6**nk);
      NNEW(inum,ITG,*nk);
      NNEW(stx,double,6*mi[0]**ne);
      NNEW(eei,double,6*mi[0]**ne);
	    
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      *iout=2;
      *icmd=3;
	    
      resultsstr(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		 ielorien,norien,orab,ntmat_,t0,t1,ithermal,
		 prestr,iprestr,filabl,eme,emn,een,iperturb,
		 f,fn,nactdof,iout,qa,vold,b,nodeboun,
		 ndirboun,xboun,nboun,ipompc,
		 nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		 bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
		 xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
		 ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		 xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		 ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		 nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
		 reltime,ne0,xforc,nforc,thicke,shcon,nshcon,
		 sideload,xload,xloadold,icfd,inomat,pslavsurf,pmastsurf,
		 mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		 islavsurf,ielprop,prop,energyini,energy,&kscale,
		 &nener,orname,&network,neapar,nebpar,ipobody,ibody,xbody,
		 nbody);
	    
      *icmd=0;
	    
      SFREE(v);SFREE(fn);SFREE(eei);
	    

      /* if the design variables are the coordinates:
	 check for the existence of a target node set */

      /* calculating the objective function */

      if(*icoordinate==1){
	nodeset=0;
	for(i=0;i<*nset;i++){
	  if(strcmp1(&objectset[m*324+162]," ")==0) continue;
	  if(strcmp2(&objectset[m*324+162],&set[i*81],81)==0){
	    nodeset=i+1;
	    break;
	  }
	}
	FORTRAN(objective_stress,(&nodeset,istartset,iendset,
				  ialset,nk,&idesvar,&iobject,mi,g0,
				  nobject,stn,objectset,&expks));
      }

      if(*icoordinate==1){
	      
	/* first design variable w.r.t. which the second derivative
	   is taken */

	for(idesvara=0;idesvara<*ndesi;idesvara++){
	      
	  NNEW(b,double,neq[1]);
		  
	  /* copying the RHS from field df */
		    
	  for(j=jqs[idesvara]-1;j<jqs[idesvara+1]-1;j++){
	    b[irows[j]-1]=df[j];
	  }
		  
	  /* solve the system */
		    
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(&b[neq[1]],&neq[1]);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(&b[neq[1]],token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(&b[neq[1]],&neq[1]);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(&b[neq[1]],&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	      
	  /* second design variable w.r.t. which the second derivative
	     is taken */

	  idesvarb=-1;num_cpusd=0;idesvarbref=0;

	  do{
	    idesvarb++;
	    num_cpusd++;
	    RENEW(b,double,(num_cpusd+1)*neq[1]);
		  
	    /* copying the RHS from field df */
			
	    for(j=jqs[idesvarb]-1;j<jqs[idesvarb+1]-1;j++){
	      b[num_cpusd*neq[1]+irows[j]-1]=df[j];
	    }
			
	    if((idesvarb==*ndesi-1)||(num_cpusd==num_cpus)){
		  
	      /* solve the system */
			    
	      if(*isolver==0){
#ifdef SPOOLES
		for(k=1;k<=num_cpusd;k++){
		  spooles_solve(&b[k*neq[1]],&neq[1]);
		}
#endif
	      }
	      else if(*isolver==4){
#ifdef SGI
		for(k=1;k<=num_cpusd;k++){
		  sgi_solve(&b[k*neq[1]],token);
		}
#endif
	      }
	      else if(*isolver==5){
#ifdef TAUCS
		for(k=1;k<=num_cpusd;k++){
		  tau_solve(&b[k*neq[1]],&neq[1]);
		}
#endif
	      }
	      else if(*isolver==7){
#ifdef PARDISO
		pardiso_solve(&b[neq[1]],&neq[1],&symmetryflag,&inputformat,&num_cpusd);
#endif
	      }

	      printf(" Using up to %" ITGFORMAT " cpu(s) for the 2nd order\n stress sensitivity.\n\n", num_cpusd);
			    
	      co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;
	      ne1=ne;stn1=stn;elcon1=elcon;nelcon1=nelcon;
	      rhcon1=rhcon;nrhcon1=nrhcon;alcon1=alcon;
	      nalcon1=nalcon;alzero1=alzero;ielmat1=ielmat;
	      ielorien1=ielorien;norien1=norien;orab1=orab;
	      ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;
	      prestr1=prestr;iprestr1=iprestr;filabl1=filabl;
	      emn1=emn;een1=een;iperturb1=iperturb;f1=f;
	      nactdof1=nactdof;vold1=vold;nodeboun1=nodeboun;
	      ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;
	      ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;
	      labmpc1=labmpc;nmpc1=nmpc;nmethod1=nmethod;
	      cam1=cam;neq1=neq;veold1=veold;accold1=accold;
	      bet1=bet;gam1=gam;dtime1=dtime;time1=time;
	      ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
	      plkcon1=plkcon;nplkcon1=nplkcon;
	      xstateini1=xstateini;xstate1=xstate;npmat1_=npmat_;
	      epn1=epn;matname1=matname;mi1=mi;ielas1=ielas;
	      ncmat1_=ncmat_;nstate1_=nstate_;stiini1=stiini;
	      vini1=vini;ikboun1=ikboun;ilboun1=ilboun;
	      enern1=enern;emeini1=emeini;xstaten1=xstaten;
	      enerini1=enerini;cocon1=cocon;ncocon1=ncocon;
	      set1=set;nset1=nset;istartset1=istartset;
	      iendset1=iendset;ialset1=ialset;nprint1=nprint;
	      prlab1=prlab;prset1=prset;qfx1=qfx;qfn1=qfn;
	      trab1=trab;inotr1=inotr;ntrans1=ntrans;fmpc1=fmpc;
	      nelemload1=nelemload;nload1=nload;ikmpc1=ikmpc;
	      ilmpc1=ilmpc;istep1=istep;iinc1=iinc;
	      springarea1=springarea;reltime1=reltime;ne01=ne0;
	      xforc1=xforc;nforc1=nforc;thicke1=thicke;
	      shcon1=shcon;nshcon1=nshcon;sideload1=sideload;
	      xload1=xload;xloadold1=xloadold;icfd1=icfd;
	      inomat1=inomat;pslavsurf1=pslavsurf;
	      pmastsurf1=pmastsurf;mortar1=mortar;
	      islavact1=islavact;cdn1=cdn;islavnode1=islavnode;
	      nslavnode1=nslavnode;ntie1=ntie;clearini1=clearini;
	      islavsurf1=islavsurf;ielprop1=ielprop;prop1=prop;
	      energyini1=energyini;energy1=energy;kscale1=kscale;
	      orname1=orname;network1=network;nestart1=nestart;
	      neend1=neend;jqs1=jqs;irows1=irows;
	      nodedesi1=nodedesi;xdesi1=xdesi;ndesi1=ndesi;
	      iobject1=iobject;nobject1=nobject;
	      objectset1=objectset;g01=g0;dgdxdy1=dgdxdy;
	      nasym1=nasym;isolver1=isolver;distmin1=distmin;
	      nodeset1=nodeset;b1=b;idesvara1=idesvara;
	      idesvarbref1=idesvarbref;ipobody1=ipobody;ibody1=ibody;
	      xbody1=xbody;nbody1=nbody;
		
	      NNEW(ithread,ITG,num_cpusd);
			    
	      /* Total difference of the mass */
	      /* create threads and wait */
			    
	      for(i=0;i<num_cpusd;i++)  {
		ithread[i]=i;
		pthread_create(&tid[i], NULL, (void *)stress_sen_2ndordermt, (void *)&ithread[i]);
	      }
			    
	      for(i=0;i<num_cpusd;i++)  pthread_join(tid[i], NULL);
			    
	      SFREE(ithread);

	      if(idesvarb==*ndesi-1){
		break;
	      }else{
		num_cpusd=0;
		idesvarbref=idesvarb;
	      }
	    }
	  }while(1);
	}
      }

      /* reactivating all elements */

      for(i=0;i<*ne;i++){
	if(ipkon[i]<-1) ipkon[i]=-2-ipkon[i];
      }

      SFREE(inum);SFREE(stn);SFREE(filabl);
      SFREE(neapar);SFREE(nebpar);
    }    
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
    pardiso_cleanup(&neq[1],&symmetryflag,&inputformat);
#endif
  }
    
  return;
    
}

/* -----------------------------------------------------------*/
/* subroutine for multithreading of the stress sensitivity    */
/* -----------------------------------------------------------*/

void *stress_sen_dx_robustmt(ITG *i){
  
  ITG idesvar,node,nea,neb,naneigh,nbneigh,neaneigh,nebneigh,j,
    node1,node2,nelem;

  /* next design variable to tread (FORTRAN-notation) */
  
  idesvar=nodedesired1[*ipos1+(*i)];
  node=nodedesi1[idesvar-1];

  nea=istartdesi1[idesvar-1];
  neb=istartdesi1[idesvar]-1;
  naneigh=istartnneigh1[node-1];
  nbneigh=istartnneigh1[node]-1;
  neaneigh=istarteneigh1[node-1];
  nebneigh=istarteneigh1[node]-1;

  memcpy(&conew1[3**nk1**i],&co1[0],sizeof(double)*3**nk1);
  memcpy(&dstx1[6*mi1[0]**ne1**i],&stx1[0],sizeof(double)*6*mi1[0]**ne1);
  memcpy(&dstn1[6**nk1**i],&stn1[0],sizeof(double)*6**nk1);

  /* pertubation of the coordinates of the design variables */ 

  for(j=0;j<3;j++){    
    conew1[(node-1)*3+j+3**nk1**i]=co1[(node-1)*3+j]+xdesi1[(idesvar-1)*3+j];
  }
  
  /* perturbation of the coordinates of the neighboring nodes of the design variables
     in case of an axisymmetric or plain strain (stress? shell?) model */
  
  nelem=ialdesi1[nea-1]-1;
  if((strcmp1(&lakon1[nelem*8+6],"A")==0)||(strcmp1(&lakon1[nelem*8+6],"E")==0)){ 
    node1=iponod2dto3d1[2*(node-1)];
    node2=iponod2dto3d1[2*(node-1)+1];
  
    for(j=0;j<3;j++){ 
      conew1[(node1-1)*3+j+3**nk1**i]=co1[(node1-1)*3+j]+xdesi1[(idesvar-1)*3+j];
      conew1[(node2-1)*3+j+3**nk1**i]=co1[(node2-1)*3+j]+xdesi1[(idesvar-1)*3+j];
    }
  } 

  stress_sen_dx(&conew1[3**nk1**i],nk1,kon1,ipkon1,lakon1,ne1,
		&dstn1[6**nk1**i],
		elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,ielmat1,
		ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
		iprestr1,filabl1,iperturb1,vold1,nmethod1,dtime1,time1,
		ttime1,plicon1,nplicon1,plkcon1,nplkcon1,xstateini1,xstate1,
		npmat1_,matname1,mi1,ielas1,ncmat1_,nstate1_,stiini1,vini1,
		emeini1,enerini1,istep1,iinc1,springarea1,reltime1,ne01,
		thicke1,pslavsurf1,pmastsurf1,mortar1,clearini1,ielprop1,
		prop1,&kscale1,&iobject1,objectset1,g01,dgdx1,
		&nea,&neb,nasym1,distmin1,&idesvar,&dstx1[6*mi1[0]**ne1**i],ialdesi1,
		ialeneigh1,&neaneigh,&nebneigh,ialnneigh1,&naneigh,&nbneigh,
		stn1,expks1,ndesi1);
  
  return NULL;
}

/* -----------------------------------------------------------*/
/* subroutine for multithreading of the stress sensitivity    */
/* -----------------------------------------------------------*/

void *stress_sen_dv_robustmt(ITG *i){
  
  ITG idof,node,nea,neb,naneigh,nbneigh,neaneigh,nebneigh,
    inode,idir,node1,node2,nelem;;

  idof=nactdofred1[*ipos1+(*i)];					  
  inode=nactdofinv1[idof];
  idir=inode-(*mt1)*(inode/(*mt1));
  node=inode/(*mt1)+1;
  
  nea=istartnk1[node-1];
  neb=istartnk1[node]-1;
  naneigh=istartnneigh1[node-1];
  nbneigh=istartnneigh1[node]-1;
  neaneigh=istarteneigh1[node-1];
  nebneigh=istarteneigh1[node]-1;
  
  memcpy(&dv1[*mt1**nk1**i],&vold1[0],sizeof(double)**mt1**nk1);
  memcpy(&dstx1[6*mi1[0]**ne1**i],&stx1[0],sizeof(double)*6*mi1[0]**ne1);
  memcpy(&dstn1[6**nk1**i],&stn1[0],sizeof(double)*6**nk1);

  /* perturbation of the coordinates of the displacements */ 

  dv1[(node-1)**mt1+idir+*mt1**nk1**i]+=dispmin1;	       

  /* perturbation of the coordinates of the neighboring nodes of the design variables
     in case of an axisymmetric or plain strain (plane stress? shell?) model */
  
  nelem=ialnk1[nea-1]-1;
  if((strcmp1(&lakon1[nelem*8+6],"A")==0)||(strcmp1(&lakon1[nelem*8+6],"E")==0)){ 
    node1=iponod2dto3d1[2*(node-1)];
    node2=iponod2dto3d1[2*(node-1)+1];
  
    dv1[(node1-1)**mt1+idir+*mt1**nk1**i]+=dispmin1;
    dv1[(node2-1)**mt1+idir+*mt1**nk1**i]+=dispmin1;
	
  } 
    
  stress_sen_dv(co1,nk1,kon1,ipkon1,lakon1,ne1,&dstn1[6**nk1**i],
		elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,ielmat1,
		ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
		iprestr1,filabl1,iperturb1,&dv1[*mt1**nk1**i],nmethod1,dtime1,
		time1,ttime1,plicon1,nplicon1,plkcon1,nplkcon1,xstateini1,
		xstate1,npmat1_,matname1,mi1,ielas1,ncmat1_,nstate1_,
		stiini1,vini1,emeini1,enerini1,istep1,iinc1,springarea1,
		reltime1,ne01,thicke1,pslavsurf1,pmastsurf1,mortar1,
		clearini1,ielprop1,prop1,&kscale1,&iobject1,g01,&nea,&neb,
		nasym1,distmin1,&dstx1[6*mi1[0]**ne1**i],ialnk1,dgdu1,ialeneigh1,
		&neaneigh,&nebneigh,ialnneigh1,&naneigh,&nbneigh,stn1,expks1,
		objectset1,&idof,&node,&idir,vold1,&dispmin1);
 
  return NULL;
}

/* -----------------------------------------------------------*/
/* subroutine for multithreading of the stress sensitivity    */
/* -----------------------------------------------------------*/

void *stress_sen_2ndordermt(ITG *i){
  
  ITG idesvarb,indexb;

  /* next design variable to tread (FORTRAN-notation) */
  
  idesvarb=idesvarbref1+(*i);
  indexb=(*i+1)*neq1[1];
  
  stress_sen_2ndorder(co1,nk1,kon1,ipkon1,lakon1,ne1,stn1,
		      elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
		      ielmat1,
		      ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,
		      prestr1,iprestr1,filabl1,emn1,een1,iperturb1,
		      f1,nactdof1,vold1,nodeboun1,
		      ndirboun1,xboun1,nboun1,ipompc1,
		      nodempc1,coefmpc1,labmpc1,nmpc1,nmethod1,cam1,neq1,
		      veold1,accold1,
		      bet1,gam1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,
		      nplkcon1,
		      xstateini1,xstate1,npmat1_,epn1,matname1,mi1,ielas1,
		      ncmat1_,nstate1_,stiini1,vini1,ikboun1,ilboun1,enern1,
		      emeini1,
		      xstaten1,enerini1,cocon1,ncocon1,set1,nset1,istartset1,
		      iendset1,
		      ialset1,nprint1,prlab1,prset1,qfx1,qfn1,trab1,inotr1,
		      ntrans1,fmpc1,
		      nelemload1,nload1,ikmpc1,ilmpc1,istep1,iinc1,springarea1,
		      reltime1,ne01,xforc1,nforc1,thicke1,shcon1,nshcon1,
		      sideload1,xload1,xloadold1,icfd1,inomat1,pslavsurf1,
		      pmastsurf1,
		      mortar1,islavact1,cdn1,islavnode1,nslavnode1,ntie1,
		      clearini1,
		      islavsurf1,ielprop1,prop1,energyini1,energy1,&kscale1,
		      orname1,&network1,&nestart1,&neend1,jqs1,irows1,nodedesi1,
		      xdesi1,ndesi1,&iobject1,nobject1,objectset1,g01,dgdx1,
		      &idesvara1,&idesvarb,nasym1,isolver1,distmin1,&nodeset1,
		      &b1[0],&b1[indexb],dgdxdy1,ipobody1,ibody1,xbody1,nbody1);
  
  return NULL;
}
